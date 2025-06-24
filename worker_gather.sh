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
    echo "Found report file: ${REPORT_FILE}"
    
    # Check if harmonization statistics already exist in the report
    if grep -q "=== HARMONIZATION STATISTICS ===" "${REPORT_FILE}"; then
        echo "Harmonization statistics already exist in report. Skipping update."
    else
        echo "Appending harmonization statistics to: ${REPORT_FILE}"
        
        # Debug: Check what log files exist
        echo "Checking for log files in gwas_log/..."
        ls -la gwas_log/scatter_*_*.out 2>/dev/null || echo "No scatter log files found"
        
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
        chromosomes_processed=0
        
        echo "Chromosome-by-chromosome breakdown:" >> "${REPORT_FILE}"
        echo "" >> "${REPORT_FILE}"
        
        # Process scatter job logs to extract harmonization statistics
        for log_file in gwas_log/scatter_*_*.out; do
            if [ -f "${log_file}" ]; then
                echo "Processing log file: ${log_file}"
                
                # Extract chromosome number from log file
                chrom=$(grep "Assigned to process chromosome:" "${log_file}" | awk '{print $5}' | head -1)
                
                if [ -n "${chrom}" ]; then
                    echo "Found chromosome: ${chrom}"
                    
                    # Extract basic statistics from the log file
                    original_count=$(grep "Original variants:" "${log_file}" | awk '{print $3}' | tail -1)
                    processed_count=$(grep "Successfully processed:" "${log_file}" | awk '{print $3}' | tail -1)
                    discarded_count=$(grep "Discarded.*unresolvable" "${log_file}" | awk '{print $3}' | tail -1)
                    discard_rate=$(grep "Discard rate:" "${log_file}" | awk '{print $3}' | tail -1)
                    
                    # Extract fixref statistics if available (look for the NS lines)
                    ref_match=$(grep "NS.*ref match" "${log_file}" | awk '{print $3}' | tail -1)
                    swapped=$(grep "NS.*swapped" "${log_file}" | awk '{print $3}' | tail -1)
                    flipped=$(grep "NS.*flipped" "${log_file}" | awk '{print $3}' | tail -1)
                    unresolved=$(grep "NS.*unresolved" "${log_file}" | awk '{print $3}' | tail -1)
                    
                    echo "Debug: original=${original_count}, processed=${processed_count}, ref_match=${ref_match}, swapped=${swapped}"
                    
                    # Only add to report if we found data
                    if [ -n "${original_count}" ] && [ "${original_count}" -gt 0 ] 2>/dev/null; then
                        chromosomes_processed=$((chromosomes_processed + 1))
                        
                        cat >> "${REPORT_FILE}" << EOF
Chromosome ${chrom}:
  Original variants: ${original_count:-0}
  Successfully processed: ${processed_count:-0}
  Discarded: ${discarded_count:-0} (${discard_rate:-N/A})
EOF
                        
                        if [ -n "${ref_match}" ] && [ "${ref_match}" -gt 0 ] 2>/dev/null; then
                            cat >> "${REPORT_FILE}" << EOF
  Reference match: ${ref_match:-0}
  REF/ALT swapped: ${swapped:-0}
  Strand flipped: ${flipped:-0}
  Unresolved: ${unresolved:-0}
EOF
                        fi
                        echo "" >> "${REPORT_FILE}"
                        
                        # Add to totals (with safety checks)
                        total_original=$((total_original + ${original_count:-0}))
                        total_processed=$((total_processed + ${processed_count:-0}))
                        total_discarded=$((total_discarded + ${discarded_count:-0}))
                        if [ -n "${ref_match}" ] && [ "${ref_match}" -gt 0 ] 2>/dev/null; then
                            total_ref_match=$((total_ref_match + ${ref_match:-0}))
                            total_swapped=$((total_swapped + ${swapped:-0}))
                            total_flipped=$((total_flipped + ${flipped:-0}))
                            total_unresolved=$((total_unresolved + ${unresolved:-0}))
                        fi
                    else
                        echo "No valid data found for chromosome ${chrom}"
                    fi
                else
                    echo "No chromosome found in log file: ${log_file}"
                fi
            fi
        done
        
        echo "Processed ${chromosomes_processed} chromosomes"
        echo "Total original variants found: ${total_original}"
        
        # Only add summary if we found data
        if [ ${total_original} -gt 0 ]; then
            # Calculate overall percentages
            processed_percent=$(echo "scale=1; ${total_processed} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
            discarded_percent=$(echo "scale=1; ${total_discarded} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
            
            if [ ${total_ref_match} -gt 0 ]; then
                ref_match_percent=$(echo "scale=1; ${total_ref_match} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
                swapped_percent=$(echo "scale=1; ${total_swapped} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
                flipped_percent=$(echo "scale=1; ${total_flipped} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
                unresolved_percent=$(echo "scale=1; ${total_unresolved} * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")
            else
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
Harmonization success rate: $(echo "scale=1; (${total_ref_match} + ${total_swapped} + ${total_flipped}) * 100 / ${total_original}" | bc -l 2>/dev/null || echo "N/A")%
Data quality: $([ $(echo "${discarded_percent} < 20" | bc -l 2>/dev/null || echo "0") -eq 1 ] && echo "Good" || echo "Moderate") (${discarded_percent}% variants discarded)

Final VCF contains ${total_processed} harmonized variants from ${total_original} original variants.
Report completed: $(date)
EOF
        else
            cat >> "${REPORT_FILE}" << EOF
=== OVERALL SUMMARY ===
No harmonization statistics found in log files.
This may indicate that the scatter jobs have not completed yet or log files are missing.
Check the gwas_log/ directory for scatter job outputs.

Report generated: $(date)
EOF
        fi
        
        echo "Harmonization statistics added to report: ${REPORT_FILE}"
    fi
else
    echo "Warning: Could not find harmonization report file to update"
fi

echo "Cleaning up temporary files..."
rm -r ${TEMP_DIR}

echo "--- Gather Job Finished Successfully ---"