# VCF Harmonization SLURM Pipeline

This pipeline harmonizes the reference alleles in a VCF file against a standard reference genome (GRCh37 or GRCh38). It is designed to run on a High-Performance Computing (HPC) cluster using the SLURM workload manager.

The process uses a scatter-gather approach:
1.  **Scatter:** The VCF is split by chromosome, and each chromosome is processed in parallel as a separate job.
2.  **Gather:** Once all chromosome jobs complete successfully, a final job merges the results into a single, harmonized VCF file.

**DISCLAIMER: This software is provided "AS IS", without warranty of any kind, express or implied. The authors and distributors of this software assume no responsibility for any loss or damage resulting from its use.**

---

### Dependencies

You must have the following software available in your cluster environment, typically loaded via a module system:

*   **SLURM:** For job submission (`sbatch`).
*   **bash:** The scripting language.
*   **BCFtools:** For VCF manipulation (`=1.14`). *Tested with bcftools 1.14.*
*   **SAMtools:** For FASTA indexing. *Tested with version 1.18 using htslib 1.18.*
*   **tabix:** For indexing compressed VCF files. *Tested with version 0.2.5 (r1005).*
*   **wget** or **curl:** For downloading reference genomes.

---

### Setup

1.  Place the three scripts (`run_vcfharmonizer.sh`, `worker_harmonize.sh`, `worker_gather.sh`) in the same directory.
2.  Make the main script executable: `chmod +x run_vcfharmonizer.sh`

---

### SLURM Configuration

The main script includes the following SLURM resource allocation settings:

*   **Job name:** `harmonization_pipeline`
*   **Output log:** `gwas_log/harmonization_pipeline-%j.out` (where %j is the job ID)
*   **Error log:** `gwas_log/harmonization_pipeline-%j.err`
*   **Time limit:** 4 hours
*   **Memory:** 32GB
*   **Tasks:** 1
*   **CPUs per task:** 4

**Note:** Ensure the `gwas_log/` directory exists before running the pipeline, or modify the output/error paths as needed for your environment.

---

### Usage

The pipeline is launched using the `run_vcfharmonizer.sh` script.