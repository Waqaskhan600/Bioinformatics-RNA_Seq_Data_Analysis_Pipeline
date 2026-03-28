#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

SECONDS=0

# ==============================================================================
# 1. Variable Definitions
# ==============================================================================
PROJECT_DIR="/media/user/New_Volume/RNA_pipeline_update"
DATA_DIR="${PROJECT_DIR}/Reads"
OUTPUT_DIR="${PROJECT_DIR}/output"

REFERENCE_INDEX="${PROJECT_DIR}/references/index/chr22.fa" # HISAT2 Index base
REFERENCE_GTF="${PROJECT_DIR}/references/gtf/chr22.gtf"    # Annotation GTF

TRIMMOMATIC_JAR="${PROJECT_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"

SAMPLE="HBR_1" # Sample prefix for processing
THREADS=6

# Dynamic sample output directory
SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"

# ==============================================================================
# 2. Folder Structure Verification & Creation
# ==============================================================================
echo "==> Setting up output directories in ${SAMPLE_OUT}..."
mkdir -p "${SAMPLE_OUT}/qc"
mkdir -p "${SAMPLE_OUT}/trimmed/qc"
mkdir -p "${SAMPLE_OUT}/alignments"
mkdir -p "${SAMPLE_OUT}/quant"
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
fastqc "${DATA_DIR}/${SAMPLE}_R1.fq" "${DATA_DIR}/${SAMPLE}_R2.fq" -o "${SAMPLE_OUT}/qc" -t ${THREADS}

# Run trimmomatic to trim reads with poor quality
echo "  -> Running Trimmomatic..."
java -jar "$TRIMMOMATIC_JAR" PE -threads ${THREADS} \
    "${DATA_DIR}/${SAMPLE}_R1.fq" "${DATA_DIR}/${SAMPLE}_R2.fq" \
    "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.paired.R1.fq" "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.unpaired.R1.fq" \
    "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.paired.R2.fq" "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.unpaired.R2.fq" \
    TRAILING:10 -phred33

echo "  -> Trimmomatic finished running for ${SAMPLE}!"

# Run fastqc on trimmed reads
echo "  -> Running FastQC on trimmed reads..."
fastqc "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.paired.R1.fq" "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.paired.R2.fq" -o "${SAMPLE_OUT}/trimmed/qc" -t ${THREADS}
echo "  -> Fastqc completed for ${SAMPLE}"

# Run HISAT2 for alignment
echo "  -> Running HISAT2 and formatting with samtools..."
hisat2 -q -x "$REFERENCE_INDEX" -1 "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.paired.R1.fq" -2 "${SAMPLE_OUT}/trimmed/${SAMPLE}.trimmed.paired.R2.fq" -p ${THREADS} | \
    samtools sort -@ ${THREADS} -o "${SAMPLE_OUT}/alignments/${SAMPLE}.bam"
echo "  -> HISAT2 finished running for ${SAMPLE}!"


echo "--------------------------------------------------"
echo "Aggregating Quality Control Reports"
echo "--------------------------------------------------"

# Run MultiQC
echo "Running MultiQC on initial FastQC..."
multiqc "${SAMPLE_OUT}/qc" -o "${SAMPLE_OUT}/qc" || echo "Warning: multiqc failed"

# Run MultiQC on Trimmed FastQC reports
echo "Running MultiQC on trimmed FastQC..."
multiqc "${SAMPLE_OUT}/trimmed/qc" -o "${SAMPLE_OUT}/trimmed/qc" || echo "Warning: multiqc failed"

echo "--------------------------------------------------"
echo "Running Quantification"
echo "--------------------------------------------------"

# Run featureCounts for quantification on the BAM file
echo "Running featureCounts..."
featureCounts -p -a "$REFERENCE_GTF" -T ${THREADS} -o "${SAMPLE_OUT}/quant/${SAMPLE}_output_data.txt" "${SAMPLE_OUT}/alignments/${SAMPLE}.bam"

# clean feature matrix
echo "Cleaning feature matrix..."
cut -f1,7- -s "${SAMPLE_OUT}/quant/${SAMPLE}_output_data.txt" > "${SAMPLE_OUT}/quant/${SAMPLE}_counts_data.txt"

echo "featureCounts finished running!"
echo "Pipeline execution completed successfully!"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."