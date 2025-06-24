#!/bin/bash

# This script runs only after all array jobs complete successfully.
set -e # Exit immediately if a command fails

echo "--- Gather Job Started ---"
echo "All worker jobs completed. Merging results."

# --- CONFIGURATION ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <final_output_vcf> <temp_dir> <chrom_style>" >&2
    exit 1
fi
FINAL_OUTPUT_VCF=$1
TEMP_DIR=$2
CHROM_STYLE=$3 # 'numeric' or 'chr'

ml load bcftools samtools

# --- MERGE ---
# Create a sorted list of the temporary VCF files to concatenate
ls ${TEMP_DIR}/*.harmonized.vcf.gz | sort -V > ${TEMP_DIR}/file_list.txt
echo "Concatenating the following files:"
cat ${TEMP_DIR}/file_list.txt

# Use a temporary name for the concatenated file
CONCAT_VCF="${TEMP_DIR}/concatenated.vcf.gz"
bcftools concat --threads ${SLURM_CPUS_PER_TASK:-1} -f ${TEMP_DIR}/file_list.txt -a -O z -o ${CONCAT_VCF}
echo "Concatenation complete."

# --- APPLY CHROMOSOME STYLE ---
if [ "${CHROM_STYLE}" == "chr" ]; then
    echo "Applying 'chr' prefix to chromosome names..."
    # Create a mapping file for bcftools
    bcftools view -h ${CONCAT_VCF} | grep "^##contig" | sed -E 's/.*ID=([^,>]+).*/\1\tchr\1/' > ${TEMP_DIR}/chr_map.txt
    
    bcftools annotate --rename-chrs ${TEMP_DIR}/chr_map.txt ${CONCAT_VCF} -O z -o ${FINAL_OUTPUT_VCF}
    echo "'chr' prefix applied."
else
    echo "Keeping 'numeric' chromosome style. No changes needed."
    mv ${CONCAT_VCF} ${FINAL_OUTPUT_VCF}
fi

# --- INDEX and CLEANUP ---
echo "Indexing final file..."
tabix -p vcf ${FINAL_OUTPUT_VCF}

echo "Final file created and indexed at: ${FINAL_OUTPUT_VCF}"

# --- UPDATE HARMONIZATION REPORT ---
echo "Updating harmonization report with detailed statistics..."

# Find the report file
INPUT_BASENAME=$(basename "${FINAL_OUTPUT_VCF}")
INPUT_BASENAME_NO_EXT="${INPUT_BASENAME%.vcf.gz}"
INPUT_BASENAME_NO_EXT="${INPUT_BASENAME_NO_EXT%_harm.allchr}"
INPUT_BASENAME_NO_EXT="${INPUT_BASENAME_NO_EXT%.allchr_test*}"
INPUT_BASENAME_NO_EXT="${INPUT_BASENAME_NO_EXT%.allchr}"
REPORT_FILE="${INPUT_BASENAME_NO_EXT}.harm.report.txt"

# Look for the report file in current directory if not found with exact name
if [ ! -f "${REPORT_FILE}" ]; then
    REPORT_FILE=$(find . -maxdepth 1 -name "*.harm.report.txt" 2>/dev/null | head -1)
fi

if [ -f "${REPORT_FILE}" ]; then
    echo "Appending harmonization statistics to: ${REPORT_FILE}"
    
    # Add harmonization statistics section
    cat >> "${REPORT_FILE}" << 'EOF'

=== HARMONIZATION STATISTICS ===
EOF
    
    # Initialize totals
    total_original=0
    total_processed=0
    total_discarded=0
    total_ref_match=0
    total_swapped=0
    total_flipped=0
    total_unresolved=0
    
    echo "Chromosome-by-chromosome breakdown:" >> "${REPORT_FILE}"
    echo "" >> "${REPORT_FILE}"
    
    # Process scatter job logs to extract harmonization statistics
    for log_file in gwas_log/scatter_*_*.out; do
        if [ -f "${log_file}" ]; then
            # Extract chromosome number from log file
            chrom=$(grep "Assigned to process chromosome:" "${log_file}" | awk '{print $5}' | head -1)
            
            if [ -n "${chrom}" ]; then
                # Extract statistics from the log file
                original_count=$(grep "Original variants:" "${log_file}" | awk '{print $3}' | tail -1)
                processed_count=$(grep "Successfully processed:" "${log_file}" | awk '{print $3}' | tail -1)
                discarded_count=$(grep "Discarded (unresolvable):" "${log_file}" | awk '{print $3}' | tail -1)
                discard_rate=$(grep "Discard rate:" "${log_file}" | awk '{print $3}' | tail -1)
                
                # Extract fixref statistics if available
                ref_match=$(grep "^NS.*ref match" "${log_file}" | awk '{print $3}' | tail -1)
                swapped=$(grep "^NS.*swapped" "${log_file}" | awk '{print $3}' | tail -1)
                flipped=$(grep "^NS.*flipped" "${log_file}" | awk '{print $3}' | tail -1)
                unresolved=$(grep "^NS.*unresolved" "${log_file}" | awk '{print $3}' | tail -1)
                
                # Only add to report if we found data
                if [ -n "${original_count}" ] && [ "${original_count}" != "0" ]; then
                    cat >> "${REPORT_FILE}" << EOF
Chromosome ${chrom}:
  Original variants: ${original_count:-0}
  Successfully processed: ${processed_count:-0}
  Discarded: ${discarded_count:-0} (${discard_rate:-N/A})
EOF
                    
                    if [ -n "${ref_match}" ] && [ "${ref_match}" != "0" ]; then
                        cat >> "${REPORT_FILE}" << EOF
  Reference match: ${ref_match:-0}
  REF/ALT swapped: ${swapped:-0}
  Strand flipped: ${flipped:-0}
  Unresolved: ${unresolved:-0}
EOF
                    fi
                    echo "" >> "${REPORT_FILE}"
                    
                    # Add to totals
                    total_original=$((total_original + ${original_count:-0}))
                    total_processed=$((total_processed + ${processed_count:-0}))
                    total_discarded=$((total_discarded + ${discarded_count:-0}))
                    total_ref_match=$((total_ref_match + ${ref_match:-0}))
                    total_swapped=$((total_swapped + ${swapped:-0}))
                    total_flipped=$((total_flipped + ${flipped:-0}))
                    total_unresolved=$((total_unresolved + ${unresolved:-0}))
                fi
            fi
        fi
    done
    
    # Calculate overall percentages
    if [ ${total_original} -gt 0 ]; then
        processed_percent=$(echo "scale=1; ${total_processed} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
        discarded_percent=$(echo "scale=1; ${total_discarded} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
        ref_match_percent=$(echo "scale=1; ${total_ref_match} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
        swapped_percent=$(echo "scale=1; ${total_swapped} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
        flipped_percent=$(echo "scale=1; ${total_flipped} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
        unresolved_percent=$(echo "scale=1; ${total_unresolved} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
    else
        processed_percent="N/A"
        discarded_percent="N/A"
        ref_match_percent="N/A"
        swapped_percent="N/A"
        flipped_percent="N/A"
        unresolved_percent="N/A"
    fi
    
    # Add overall summary
    cat >> "${REPORT_FILE}" << EOF
=== OVERALL SUMMARY ===
Total original variants: ${total_original}
Total successfully processed: ${total_processed} (${processed_percent}%)
Total discarded: ${total_discarded} (${discarded_percent}%)

Detailed breakdown:
  Reference match (no changes needed): ${total_ref_match} (${ref_match_percent}%)
  REF/ALT swapped (corrected): ${total_swapped} (${swapped_percent}%)
  Strand flipped (corrected): ${total_flipped} (${flipped_percent}%)
  Unresolved (discarded): ${total_unresolved} (${unresolved_percent}%)

=== HARMONIZATION QUALITY ===
EOF

    # Add quality assessment
    if [ ${total_original} -gt 0 ]; then
        success_rate=$(echo "scale=1; (${total_ref_match} + ${total_swapped} + ${total_flipped}) * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
        quality_rating="Good"
        if [ $(echo "${discarded_percent} > 20" | bc -l 2>/dev/null) -eq 1 ]; then
            quality_rating="Moderate"
        fi
        if [ $(echo "${discarded_percent} > 40" | bc -l 2>/dev/null) -eq 1 ]; then
            quality_rating="Poor"
        fi
        
        cat >> "${REPORT_FILE}" << EOF
Harmonization success rate: ${success_rate}%
Data quality: ${quality_rating} (${discarded_percent}% variants discarded)

Interpretation:
- Reference match: Variants that already matched the reference genome
- REF/ALT swapped: Variants where REF and ALT alleles were swapped (now corrected)
- Strand flipped: Variants on opposite strand (now corrected)  
- Unresolved: Variants that could not be automatically harmonized (removed)

Final VCF contains ${total_processed} harmonized variants from ${total_original} original variants.
Report completed: $(date)
EOF
    fi
    
    echo "Harmonization statistics added to report: ${REPORT_FILE}"
else
    echo "Warning: Could not find harmonization report file to update"
fi

echo "Cleaning up temporary files..."
rm -r ${TEMP_DIR}

echo "--- Gather Job Finished Successfully ---"