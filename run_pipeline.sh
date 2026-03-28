#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

SECONDS=0

# ==============================================================================
# 1. Variable Definitions
# ==============================================================================
PROJECT_DIR="/media/user/New_Volume/RNA_pipeline_update"
DATA_DIR="${PROJECT_DIR}/Reads"
RESULTS_DIR="${PROJECT_DIR}/results"

REFERENCE_INDEX="${PROJECT_DIR}/references/index/chr22.fa" # HISAT2 Index base
REFERENCE_GTF="${PROJECT_DIR}/references/gtf/chr22.gtf"    # Annotation GTF

TRIMMOMATIC_JAR="${PROJECT_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"

SAMPLE="demo_1" # Sample prefix for processing
THREADS=6

# ==============================================================================
# 2. Folder Structure Verification & Creation
# ==============================================================================
echo "==> Setting up output directories in ${RESULTS_DIR}..."
mkdir -p "${RESULTS_DIR}/qc/fastqc1"
mkdir -p "${RESULTS_DIR}/qc/fastqc2"
mkdir -p "${RESULTS_DIR}/qc/multiqc1"
mkdir -p "${RESULTS_DIR}/qc/multiqc2"
mkdir -p "${RESULTS_DIR}/trimmed"
mkdir -p "${RESULTS_DIR}/alignments"
mkdir -p "${RESULTS_DIR}/quant"
echo "Output Directories created/verified."

# ==============================================================================
# 3. Required Packages Validation
# ==============================================================================
echo "==> Validating required packages..."
REQUIRED_CMDS=("fastqc" "java" "hisat2" "samtools" "multiqc" "featureCounts")

for cmd in "${REQUIRED_CMDS[@]}"; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: Required command '$cmd' is not installed or not in PATH."
        exit 1
    fi
done
echo "All required packages are available."

# ==============================================================================
# 4. Reference Files Validation
# ==============================================================================
echo "==> Validating reference and tool files..."
if [[ ! -f "$TRIMMOMATIC_JAR" ]]; then
    echo "Error: Trimmomatic JAR not found at $TRIMMOMATIC_JAR"
    exit 1
fi

if [[ ! -f "$REFERENCE_GTF" ]]; then
    echo "Error: GTF annotation not found at $REFERENCE_GTF"
    exit 1
fi

# HISAT2 indices usually have extensions like .1.ht2 or .1.ht2l. We check if the base prefix files exist.
if ! ls ${REFERENCE_INDEX}* &> /dev/null; then
    echo "Error: HISAT2 index files missing for base prefix $REFERENCE_INDEX"
    exit 1
fi

if [[ ! -f "${DATA_DIR}/${SAMPLE}_R1.fq" ]] || [[ ! -f "${DATA_DIR}/${SAMPLE}_R2.fq" ]]; then
    echo "Error: Paired-end fastq files for sample '${SAMPLE}' not found in ${DATA_DIR}"
    exit 1
fi
echo "Reference and input files validated successfully."

# ==============================================================================
# 5. Pipeline Execution
# ==============================================================================
echo "==> Starting pipeline execution..."

echo "--------------------------------------------------"
echo "Processing sample: ${SAMPLE}"
echo "--------------------------------------------------"

# Run fastqc on original reads
echo "  -> Running FastQC on original reads..."
fastqc "${DATA_DIR}/${SAMPLE}_R1.fq" "${DATA_DIR}/${SAMPLE}_R2.fq" -o "${RESULTS_DIR}/qc/fastqc1" -t ${THREADS}

# Run trimmomatic to trim reads with poor quality
echo "  -> Running Trimmomatic..."
java -jar "$TRIMMOMATIC_JAR" PE -threads ${THREADS} \
    "${DATA_DIR}/${SAMPLE}_R1.fq" "${DATA_DIR}/${SAMPLE}_R2.fq" \
    "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.paired.R1.fq" "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.unpaired.R1.fq" \
    "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.paired.R2.fq" "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.unpaired.R2.fq" \
    TRAILING:10 -phred33

echo "  -> Trimmomatic finished running for ${SAMPLE}!"

# Run fastqc on trimmed reads
echo "  -> Running FastQC on trimmed reads..."
fastqc "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.paired.R1.fq" "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.paired.R2.fq" -o "${RESULTS_DIR}/qc/fastqc2" -t ${THREADS}
echo "  -> Fastqc completed for ${SAMPLE}"

# Run HISAT2 for alignment
echo "  -> Running HISAT2 and formatting with samtools..."
hisat2 -q -x "$REFERENCE_INDEX" -1 "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.paired.R1.fq" -2 "${RESULTS_DIR}/trimmed/${SAMPLE}.trimmed.paired.R2.fq" -p ${THREADS} | \
    samtools sort -@ ${THREADS} -o "${RESULTS_DIR}/alignments/${SAMPLE}.bam"
echo "  -> HISAT2 finished running for ${SAMPLE}!"


echo "--------------------------------------------------"
echo "Aggregating Quality Control Reports"
echo "--------------------------------------------------"

# Run MultiQC
echo "Running MultiQC on initial FastQC..."
multiqc "${RESULTS_DIR}/qc/fastqc1" -o "${RESULTS_DIR}/qc/multiqc1" || echo "Warning: multiqc failed"

# Run MultiQC on second Fastqc reports
echo "Running MultiQC on trimmed FastQC..."
multiqc "${RESULTS_DIR}/qc/fastqc2" -o "${RESULTS_DIR}/qc/multiqc2" || echo "Warning: multiqc failed"

echo "--------------------------------------------------"
echo "Running Quantification"
echo "--------------------------------------------------"

# Run featureCounts for quantification on the BAM file
echo "Running featureCounts..."
featureCounts -p -a "$REFERENCE_GTF" -T ${THREADS} -o "${RESULTS_DIR}/quant/output_data.txt" "${RESULTS_DIR}/alignments/${SAMPLE}.bam"

# clean feature matrix
echo "Cleaning feature matrix..."
cut -f1,7- -s "${RESULTS_DIR}/quant/output_data.txt" > "${RESULTS_DIR}/quant/counts_data.txt"

echo "featureCounts finished running!"
echo "Pipeline execution completed successfully!"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."