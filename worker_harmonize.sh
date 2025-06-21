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
# The norm command harmonizes against the numeric-only reference.
# The annotate command ensures the output contig is just the number (e.g., '22').
# The sort command prepares for indexing.
echo "Running harmonization and sorting pipeline..."

bcftools view -r ${CHROM_QUERY} --threads ${SLURM_CPUS_PER_TASK:-1} ${INPUT_VCF} | \
  bcftools norm -m - --check-ref f -f ${REF_FASTA} | \
  bcftools annotate --set-id +'%CHROM' | sed "s/^chr${CHROM}/${CHROM}/" | \
  bcftools sort -O z -o ${TEMP_VCF_OUT_GZ}

echo "Harmonization and sorting command finished."

# --- FINAL CHECK AND REPORTING ---
if [ -s "${TEMP_VCF_OUT_GZ}" ]; then
    echo "[SUCCESS] Output file ${TEMP_VCF_OUT_GZ} was created successfully."

    echo "--- Generating harmonization summary report ---"
    bcftools +fixref - -- -f ${REF_FASTA} < ${TEMP_VCF_OUT_GZ} > /dev/null
    echo "--- Report finished ---"

    echo "Indexing the output file..."
    tabix -p vcf ${TEMP_VCF_OUT_GZ}
    echo "Indexing complete."
else
    echo "FATAL ERROR: Output file ${TEMP_VCF_OUT_GZ} is empty or was not created." >&2
    exit 1
fi

echo "--- Worker Job for ${CHROM} Finished Successfully ---"