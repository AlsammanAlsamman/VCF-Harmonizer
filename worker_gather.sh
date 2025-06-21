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
echo "Cleaning up temporary files..."
rm -r ${TEMP_DIR}

echo "--- Gather Job Finished Successfully ---"