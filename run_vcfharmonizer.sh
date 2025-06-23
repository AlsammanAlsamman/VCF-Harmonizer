#!/bin/bash

#SBATCH --job-name=harmonization_pipeline
#SBATCH --output=gwas_log/harmonization_pipeline-%j.out
#SBATCH --error=gwas_log/harmonization_pipeline-%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# ==============================================================================
#                      VCF Harmonization Pipeline
#
# This script orchestrates a scatter-gather SLURM pipeline to harmonize a VCF
# file against a specified reference genome.
#
# It can use a pre-existing reference FASTA or automatically download, format,
# and index a standard build (GRCh37/GRCh38) if one is not provided.
# ==============================================================================

set -e
set -o pipefail

# --- Default Configuration ---
GENOME_BUILD="GRCh37"
CHROM_STYLE="numeric" # 'numeric' (1, 2) or 'chr' (chr1, chr2)
INPUT_VCF=""
FINAL_OUTPUT_VCF=""
CUSTOM_REF_FASTA=""
REF_BASE_DIR="reference_genomes"

# --- Functions ---

# Function to display usage information
usage() {
  echo "Usage: $0 -i <input.vcf.gz> -o <output.vcf.gz> [OPTIONS]"
  echo "Options:"
  echo "  -i <file>     (Required) Input VCF file (must be bgzipped and indexed)."
  echo "  -o <file>     (Required) Final output VCF file name."
  echo "  -r <file>     Path to a custom reference FASTA file. If provided, -g is ignored."
  echo "                The FASTA should use numeric chromosome names (1, 22, X) for best results."
  echo "  -g <build>    Genome build if not providing a custom reference. Either 'GRCh37' or 'GRCh38'."
  echo "                (Default: GRCh37)"
  echo "  -c <style>    Chromosome naming style for the output VCF. Either 'numeric' or 'chr'."
  echo "                (Default: numeric)"
  echo "  -h            Display this help message."
  exit 1
}

# Function to download and prepare a reference genome
setup_reference() {
    local build=$1
    local ref_dir="${REF_BASE_DIR}/${build}"
    local final_fasta_path="${ref_dir}/${build}.numeric.fasta"

    # If the final, processed FASTA already exists, we're done.
    if [ -f "${final_fasta_path}" ]; then
        echo "Found prepared reference genome: ${final_fasta_path}"
        REF_FASTA_PATH=${final_fasta_path}
        return
    fi

    echo "--- Preparing reference genome ${build} ---"
    mkdir -p "${ref_dir}"
    
    local original_dir=$(pwd)
    cd "${ref_dir}"

    local download_url=""
    local raw_fasta_gz=""

    if [ "${build}" == "GRCh37" ]; then
        download_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz"
        raw_fasta_gz="GCF_000001405.13_GRCh37_genomic.fna.gz"
    elif [ "${build}" == "GRCh38" ]; then
        download_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz"
        raw_fasta_gz="GCF_000001405.26_GRCh38_genomic.fna.gz"
    else
        echo "Error: Invalid genome build '${build}'." >&2
        cd "${original_dir}"
        exit 1
    fi

    if [ ! -f "${raw_fasta_gz}" ]; then
        echo "Downloading ${build} from NCBI..."
        wget -O "${raw_fasta_gz}" "${download_url}"
    fi

    echo "Reformatting FASTA headers to numeric style and keeping primary chromosomes..."
    if [ "${build}" == "GRCh37" ]; then
        gunzip -c "${raw_fasta_gz}" | sed -E \
            -e 's/^>NC_000001\.10.*/>1/'   -e 's/^>NC_000002\.11.*/>2/' \
            -e 's/^>NC_000003\.11.*/>3/'   -e 's/^>NC_000004\.11.*/>4/' \
            -e 's/^>NC_000005\.9.*/>5/'    -e 's/^>NC_000006\.11.*/>6/' \
            -e 's/^>NC_000007\.13.*/>7/'   -e 's/^>NC_000008\.10.*/>8/' \
            -e 's/^>NC_000009\.11.*/>9/'   -e 's/^>NC_000010\.10.*/>10/' \
            -e 's/^>NC_000011\.9.*/>11/'   -e 's/^>NC_000012\.11.*/>12/' \
            -e 's/^>NC_000013\.10.*/>13/'   -e 's/^>NC_000014\.8.*/>14/' \
            -e 's/^>NC_000015\.9.*/>15/'   -e 's/^>NC_000016\.9.*/>16/' \
            -e 's/^>NC_000017\.10.*/>17/'   -e 's/^>NC_000018\.9.*/>18/' \
            -e 's/^>NC_000019\.9.*/>19/'   -e 's/^>NC_000020\.10.*/>20/' \
            -e 's/^>NC_000021\.8.*/>21/'   -e 's/^>NC_000022\.10.*/>22/' \
            -e 's/^>NC_000023\.10.*/>X/'   -e 's/^>NC_000024\.9.*/>Y/' \
            -e 's/^>NC_012920\.1.*/>MT/' | awk '/^>/ { if (p) print p; p=0; if ($0 ~ /^>[0-9XYM]/) p=1; } p' > "${final_fasta_path}"
    else # GRCh38
        gunzip -c "${raw_fasta_gz}" | sed -E \
            -e 's/^>NC_00000([1-9])\.11.*/>\1/' \
            -e 's/^>NC_0000(1[0-9])\.11.*/>\1/' -e 's/^>NC_0000(2[0-2])\.(11|12|10|9).*/>\1/' \
            -e 's/^>NC_000023\.11.*/>X/'       -e 's/^>NC_000024\.10.*/>Y/' \
            -e 's/^>NC_012920\.1.*/>MT/' | awk '/^>/ { if (p) print p; p=0; if ($0 ~ /^>[0-9XYM]/) p=1; } p' > "${final_fasta_path}"
    fi
    
    echo "Indexing reference..."
    samtools faidx "${final_fasta_path}"

    cd "${original_dir}"
    echo "--- Reference setup complete ---"
    REF_FASTA_PATH=${final_fasta_path}
}

# --- Argument Parsing ---
while getopts "i:o:r:g:c:h" opt; do
  case ${opt} in
    i) INPUT_VCF=$OPTARG ;;
    o) FINAL_OUTPUT_VCF=$OPTARG ;;
    r) CUSTOM_REF_FASTA=$OPTARG ;;
    g) GENOME_BUILD=$OPTARG ;;
    c) CHROM_STYLE=$OPTARG ;;
    h) usage ;;
    \?) usage ;;
  esac
done

if [ -z "${INPUT_VCF}" ] || [ -z "${FINAL_OUTPUT_VCF}" ]; then
  echo "Error: Input and output files must be specified."
  usage
fi

# Generate report file name based on input VCF
INPUT_BASENAME=$(basename "${INPUT_VCF}")
INPUT_BASENAME_NO_EXT="${INPUT_BASENAME%.vcf.gz}"
INPUT_BASENAME_NO_EXT="${INPUT_BASENAME_NO_EXT%.vcf}"
REPORT_FILE="${INPUT_BASENAME_NO_EXT}.harm.report.txt"

# --- Main Logic ---
# Load modules
echo "--- Loading required modules ---"
ml load bcftools samtools

# Prepare reference genome
if [ -n "${CUSTOM_REF_FASTA}" ]; then
    echo "--- Using custom reference file: ${CUSTOM_REF_FASTA} ---"
    if [ ! -f "${CUSTOM_REF_FASTA}" ]; then
        echo "Error: Custom reference file not found at ${CUSTOM_REF_FASTA}" >&2
        exit 1
    fi
    # Ensure it's indexed
    if [ ! -f "${CUSTOM_REF_FASTA}.fai" ]; then
        echo "Reference index (.fai) not found. Creating one..."
        samtools faidx "${CUSTOM_REF_FASTA}"
        echo "Index created at ${CUSTOM_REF_FASTA}.fai"
    fi
    REF_FASTA_PATH=${CUSTOM_REF_FASTA}
else
    echo "--- No custom reference provided. Using genome build: ${GENOME_BUILD} ---"
    setup_reference "${GENOME_BUILD}"
fi

# Setup temporary directory
TEMP_DIR="harmonization_temp_$$"
mkdir -p ${TEMP_DIR}
mkdir -p gwas_log
echo "--- Temporary directory created at ${TEMP_DIR} ---"

# Extract chromosome list from input VCF (numeric only, as our default references are numeric)
echo "Extracting chromosome list..."
bcftools view -h ${INPUT_VCF} | grep "^##contig" | sed -E 's/.*ID=(chr)?([^,>]+).*/\2/' | grep -E '^[0-9]+$|^X$|^Y$' | sort -V > ${TEMP_DIR}/chr_list.txt

N_CHROMS=$(wc -l < ${TEMP_DIR}/chr_list.txt)
if [ "$N_CHROMS" -eq 0 ]; then
    echo "Error: No standard chromosomes (1-22, X, Y) found in the input VCF file." >&2
    exit 1
fi
echo "Found ${N_CHROMS} chromosomes to process."

# --- STEP 1 & 2: SUBMIT SLURM JOBS ---
echo "Submitting parallel harmonization worker jobs..."
ARRAY_JOB_ID=$(sbatch \
  --parsable \
  --job-name=Harmonize_Scatter \
  --output=gwas_log/scatter_%A_%a.out \
  --error=gwas_log/scatter_%A_%a.err \
  --array=1-${N_CHROMS} \
  --time=02:00:00 \
  --mem=8G \
  --cpus-per-task=4 \
  worker_harmonize.sh "${INPUT_VCF}" "${REF_FASTA_PATH}" "${TEMP_DIR}")

if [ $? -ne 0 ]; then
    echo "Error submitting array job. Exiting." >&2
    exit 1
fi
echo "Array job submitted with ID: ${ARRAY_JOB_ID}"

echo "Submitting the final gather job..."
sbatch \
  --job-name=Harmonize_Gather \
  --output=gwas_log/gather_%j.out \
  --error=gwas_log/gather_%j.err \
  --dependency=afterok:${ARRAY_JOB_ID} \
  --time=01:00:00 \
  --mem=16G \
  --cpus-per-task=4 \
  worker_gather.sh "${FINAL_OUTPUT_VCF}" "${TEMP_DIR}" "${CHROM_STYLE}"

echo "--- Pipeline Submitted ---"

# Generate report file
cat > "${REPORT_FILE}" << EOF
VCF Harmonization Pipeline Report
Generated: $(date)

=== INPUT PARAMETERS ===
Input VCF: ${INPUT_VCF}
Output VCF: ${FINAL_OUTPUT_VCF}
Genome Build: ${GENOME_BUILD}
Chromosome Style: ${CHROM_STYLE}
EOF

if [ -n "${CUSTOM_REF_FASTA}" ]; then
    echo "Reference: Custom file at ${REF_FASTA_PATH}"
    echo "Reference: Custom file at ${REF_FASTA_PATH}" >> "${REPORT_FILE}"
else
    echo "Reference: Build ${GENOME_BUILD}"
    echo "Reference: Build ${GENOME_BUILD}" >> "${REPORT_FILE}"
fi

cat >> "${REPORT_FILE}" << EOF

=== JOB INFORMATION ===
Scatter Job ID: ${scatter_job_id}
Gather Job ID: ${gather_job_id}
Temporary Directory: ${TEMP_DIR}

=== MONITORING ===
Monitor jobs with: squeue -u \$USER
Check logs in: gwas_log/ directory

=== OUTPUT ===
Final harmonized VCF: ${FINAL_OUTPUT_VCF}
Report file: ${REPORT_FILE}
EOF

echo "Output Chromosome Style: ${CHROM_STYLE}"
echo "Monitor jobs with 'squeue -u \$USER'"
echo "Final output will be at: ${FINAL_OUTPUT_VCF}"
echo "Report written to: ${REPORT_FILE}"