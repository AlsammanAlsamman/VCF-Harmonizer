#!/bin/bash
# ===================================================================
#           Worker Script (No changes needed)
# ===================================================================

set -e
set -o pipefail

echo "--- Worker Job Started ---"

# --- CONFIGURATION ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_vcf> <ref_fasta> <temp_dir>"
    exit 1
fi
INPUT_VCF=$1
REF_FASTA=$2
TEMP_DIR=$3

# --- SETUP ---
ml load bcftools samtools
# Note: The chr_list.txt now only contains numeric chromosomes
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${TEMP_DIR}/chr_list.txt)
TEMP_VCF_OUT_GZ="${TEMP_DIR}/${CHROM}.harmonized.vcf.gz"

# We must query the input VCF with both 'chr' and 'numeric' styles
# because we don't know its original format.
CHROM_QUERY="${CHROM},chr${CHROM}"

echo "Assigned to process chromosome: ${CHROM} (Querying as ${CHROM_QUERY})"
echo "Reference FASTA: ${REF_FASTA}"
echo "Output will be written to: ${TEMP_VCF_OUT_GZ}"

# --- HARMONIZE, SORT, and NORMALIZE CHROMOSOME NAME ---
# The view command extracts the chromosome (e.g., 'chr22' or '22').
# The fixref plugin handles strand flips and reference mismatches automatically.
# The norm command normalizes variants and splits multiallelic sites.
# The annotate command ensures the output contig is just the number (e.g., '22').
# The sort command prepares for indexing.
echo "Running strand-flip aware harmonization pipeline..."

# Step 1: Extract chromosome and fix reference alleles (handles strand flips)
echo "Step 1: Extracting chromosome ${CHROM} and fixing reference alleles..."

# First run: Get stats to see what we're dealing with
echo "Running fixref stats analysis..."
bcftools view -r ${CHROM_QUERY} --threads ${SLURM_CPUS_PER_TASK:-1} ${INPUT_VCF} | \
  bcftools +fixref -- -f ${REF_FASTA} -m stats > ${TEMP_DIR}/${CHROM}_fixref_stats.log 2>&1

# Second run: Apply fixes, discarding unresolvable variants
echo "Applying fixref corrections..."
bcftools view -r ${CHROM_QUERY} --threads ${SLURM_CPUS_PER_TASK:-1} ${INPUT_VCF} | \
  bcftools +fixref -- -f ${REF_FASTA} -m flip -d 2> ${TEMP_DIR}/${CHROM}_fixref.log | \
  bcftools view -O z -o ${TEMP_DIR}/${CHROM}_fixed.vcf.gz

# Index the fixed file
tabix -p vcf ${TEMP_DIR}/${CHROM}_fixed.vcf.gz

echo "Fixref processing completed."

# Step 2: Normalize and finalize (with more permissive settings for remaining issues)
echo "Step 2: Normalizing and finalizing chromosome ${CHROM}..."
bcftools norm -m - --check-ref w -f ${REF_FASTA} ${TEMP_DIR}/${CHROM}_fixed.vcf.gz 2> ${TEMP_DIR}/${CHROM}_norm_final.log | \
  bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' | \
  sed "s/^chr${CHROM}/${CHROM}/" | \
  bcftools sort -O z -o ${TEMP_VCF_OUT_GZ}

# Clean up intermediate file
rm -f ${TEMP_DIR}/${CHROM}_fixed.vcf.gz ${TEMP_DIR}/${CHROM}_fixed.vcf.gz.tbi

echo "Strand-flip aware harmonization completed."

# Report fixref results
echo "=== FIXREF ANALYSIS FOR CHROMOSOME ${CHROM} ==="

# Show the stats analysis first
if [ -f "${TEMP_DIR}/${CHROM}_fixref_stats.log" ] && [ -s "${TEMP_DIR}/${CHROM}_fixref_stats.log" ]; then
    echo "=== ORIGINAL VARIANT ANALYSIS ==="
    cat "${TEMP_DIR}/${CHROM}_fixref_stats.log"
    echo ""
fi

# Calculate before/after counts
if [ -f "${TEMP_DIR}/${CHROM}_fixed.vcf.gz" ]; then
    original_count=$(bcftools view -r ${CHROM_QUERY} -H ${INPUT_VCF} | wc -l)
    fixed_count=$(bcftools view -H ${TEMP_DIR}/${CHROM}_fixed.vcf.gz | wc -l)
    discarded_count=$((original_count - fixed_count))
    
    echo "=== FIXREF RESULTS SUMMARY ==="
    echo "Original variants: ${original_count}"
    echo "Successfully processed: ${fixed_count}"
    echo "Discarded (unresolvable): ${discarded_count}"
    
    if [ ${original_count} -gt 0 ]; then
        discard_percent=$(echo "scale=1; ${discarded_count} * 100 / ${original_count}" | bc -l 2>/dev/null || echo "N/A")
        echo "Discard rate: ${discard_percent}%"
    fi
    echo ""
fi

# Show any error messages
if [ -f "${TEMP_DIR}/${CHROM}_fixref.log" ] && [ -s "${TEMP_DIR}/${CHROM}_fixref.log" ]; then
    echo "=== FIXREF PROCESSING LOG ==="
    cat "${TEMP_DIR}/${CHROM}_fixref.log"
    echo ""
fi

echo "=== END FIXREF ANALYSIS ==="

# --- FINAL CHECK AND REPORTING ---
if [ -s "${TEMP_VCF_OUT_GZ}" ]; then
    echo "[SUCCESS] Output file ${TEMP_VCF_OUT_GZ} was created successfully."
    
    # Get variant count for reporting
    variant_count=$(bcftools view -H ${TEMP_VCF_OUT_GZ} | wc -l)
    echo "Final variant count for chromosome ${CHROM}: ${variant_count}"

    echo "Indexing the output file..."
    tabix -p vcf ${TEMP_VCF_OUT_GZ}
    echo "Indexing complete."
    
    # Report final normalization issues (if any)
    if [ -f "${TEMP_DIR}/${CHROM}_norm_final.log" ] && [ -s "${TEMP_DIR}/${CHROM}_norm_final.log" ]; then
        echo "=== FINAL NORMALIZATION LOG FOR CHROMOSOME ${CHROM} ==="
        cat "${TEMP_DIR}/${CHROM}_norm_final.log"
        echo "=== END FINAL LOG ==="
    fi
else
    echo "FATAL ERROR: Output file ${TEMP_VCF_OUT_GZ} is empty or was not created." >&2
    exit 1
fi

echo "--- Worker Job for ${CHROM} Finished Successfully ---"