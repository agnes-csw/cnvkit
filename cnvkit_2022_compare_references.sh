#!/bin/bash
#
# CNVkit Pipeline for 2022 Dataset (Tumor-only samples)
# Compares two reference approaches:
#   1. Tumor-based reference (built from tumor samples)
#   2. Flat reference (neutral copy number)
#

set -e

# =============================================================================
# Configuration
# =============================================================================
export FASTA="/insomnia001/depts/pmg/users/sc5006/reference/Homo_sapiens_assembly38.fasta"
export BAM_DIR="/insomnia001/depts/pmg/users/sc5006/data/final_bams"
SAMPLE_FILE="samples.txt"
export REF_DIR="ref"

# Output directories for each approach
export OUTPUT_TUMOR_REF="cnvkit_output_tumor_ref"
export OUTPUT_FLAT_REF="cnvkit_output_flat_ref"

# Apptainer/Singularity settings
SIF_IMAGE="/manitou-home/pmg/users/sc5006/data/singularity_images/cnvkit.sif"
BIND_PATHS="-B /manitou-home -B /insomnia001"
export CNVKIT="apptainer exec $BIND_PATHS $SIF_IMAGE cnvkit.py"

# Number of parallel jobs
NJOBS=8

# Target/antitarget BED files (Twist panel)
export TARGET_BED="$REF_DIR/Twist_Comprehensive_Exome_Covered_Targets_hg38.target.bed"
export ANTITARGET_BED="$REF_DIR/Twist_Comprehensive_Exome_Covered_Targets_hg38.antitarget.bed"

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$OUTPUT_TUMOR_REF" "$OUTPUT_FLAT_REF"

# Read sample IDs into array
mapfile -t SAMPLES < "$SAMPLE_FILE"
echo "Found ${#SAMPLES[@]} tumor samples"

# =============================================================================
# Step 1: Run coverage for ALL samples (shared between both approaches)
# =============================================================================
echo "=============================================="
echo "Step 1: Running coverage analysis..."
echo "=============================================="

run_coverage() {
    sample_id=$1
    bam_file="$BAM_DIR/${sample_id}_final.bam"

    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file" >&2
        return 1
    fi

    # Target coverage (store in tumor_ref output, will be used by both)
    if [[ ! -f "$OUTPUT_TUMOR_REF/${sample_id}.targetcoverage.cnn" ]]; then
        $CNVKIT coverage "$bam_file" "$TARGET_BED" \
            -o "$OUTPUT_TUMOR_REF/${sample_id}.targetcoverage.cnn"
    fi

    # Antitarget coverage
    if [[ ! -f "$OUTPUT_TUMOR_REF/${sample_id}.antitargetcoverage.cnn" ]]; then
        $CNVKIT coverage "$bam_file" "$ANTITARGET_BED" \
            -o "$OUTPUT_TUMOR_REF/${sample_id}.antitargetcoverage.cnn"
    fi

    echo "Completed coverage: $sample_id"
}
export -f run_coverage

printf '%s\n' "${SAMPLES[@]}" | parallel -j "$NJOBS" run_coverage {}

# =============================================================================
# Step 2A: Build TUMOR-BASED reference
# =============================================================================
echo "=============================================="
echo "Step 2A: Building tumor-based reference..."
echo "=============================================="

if [[ ! -f "$REF_DIR/reference_tumor_2022.cnn" ]]; then
    # Collect all coverage files
    CNN_FILES=""
    for sample_id in "${SAMPLES[@]}"; do
        CNN_FILES="$CNN_FILES $OUTPUT_TUMOR_REF/${sample_id}.targetcoverage.cnn"
        CNN_FILES="$CNN_FILES $OUTPUT_TUMOR_REF/${sample_id}.antitargetcoverage.cnn"
    done

    $CNVKIT reference $CNN_FILES \
        --fasta "$FASTA" \
        -o "$REF_DIR/reference_tumor_2022.cnn"
fi

# =============================================================================
# Step 2B: Build FLAT reference
# =============================================================================
echo "=============================================="
echo "Step 2B: Building flat reference..."
echo "=============================================="

if [[ ! -f "$REF_DIR/reference_flat_2022.cnn" ]]; then
    $CNVKIT reference \
        -o "$REF_DIR/reference_flat_2022.cnn" \
        -f "$FASTA" \
        -t "$TARGET_BED" \
        -a "$ANTITARGET_BED"
fi

# =============================================================================
# Step 3A: Fix and segment using TUMOR-BASED reference
# =============================================================================
echo "=============================================="
echo "Step 3A: Fix/segment with tumor-based reference..."
echo "=============================================="

run_fix_segment_tumor() {
    sample_id=$1

    # Fix
    if [[ ! -f "$OUTPUT_TUMOR_REF/${sample_id}.cnr" ]]; then
        $CNVKIT fix \
            "$OUTPUT_TUMOR_REF/${sample_id}.targetcoverage.cnn" \
            "$OUTPUT_TUMOR_REF/${sample_id}.antitargetcoverage.cnn" \
            "$REF_DIR/reference_tumor_2022.cnn" \
            -o "$OUTPUT_TUMOR_REF/${sample_id}.cnr"
    fi

    # Segment
    if [[ ! -f "$OUTPUT_TUMOR_REF/${sample_id}.cns" ]]; then
        $CNVKIT segment \
            "$OUTPUT_TUMOR_REF/${sample_id}.cnr" \
            -o "$OUTPUT_TUMOR_REF/${sample_id}.cns"
    fi

    # Call
    if [[ ! -f "$OUTPUT_TUMOR_REF/${sample_id}.call.cns" ]]; then
        $CNVKIT call \
            "$OUTPUT_TUMOR_REF/${sample_id}.cns" \
            -o "$OUTPUT_TUMOR_REF/${sample_id}.call.cns"
    fi

    echo "Completed tumor-ref: $sample_id"
}
export -f run_fix_segment_tumor

printf '%s\n' "${SAMPLES[@]}" | parallel -j "$NJOBS" run_fix_segment_tumor {}

# =============================================================================
# Step 3B: Fix and segment using FLAT reference
# =============================================================================
echo "=============================================="
echo "Step 3B: Fix/segment with flat reference..."
echo "=============================================="

run_fix_segment_flat() {
    sample_id=$1

    # Fix (use same coverage files, different reference)
    if [[ ! -f "$OUTPUT_FLAT_REF/${sample_id}.cnr" ]]; then
        $CNVKIT fix \
            "$OUTPUT_TUMOR_REF/${sample_id}.targetcoverage.cnn" \
            "$OUTPUT_TUMOR_REF/${sample_id}.antitargetcoverage.cnn" \
            "$REF_DIR/reference_flat_2022.cnn" \
            -o "$OUTPUT_FLAT_REF/${sample_id}.cnr"
    fi

    # Segment
    if [[ ! -f "$OUTPUT_FLAT_REF/${sample_id}.cns" ]]; then
        $CNVKIT segment \
            "$OUTPUT_FLAT_REF/${sample_id}.cnr" \
            -o "$OUTPUT_FLAT_REF/${sample_id}.cns"
    fi

    # Call
    if [[ ! -f "$OUTPUT_FLAT_REF/${sample_id}.call.cns" ]]; then
        $CNVKIT call \
            "$OUTPUT_FLAT_REF/${sample_id}.cns" \
            -o "$OUTPUT_FLAT_REF/${sample_id}.call.cns"
    fi

    echo "Completed flat-ref: $sample_id"
}
export -f run_fix_segment_flat

printf '%s\n' "${SAMPLES[@]}" | parallel -j "$NJOBS" run_fix_segment_flat {}

# =============================================================================
# Step 4: Generate comparison metrics
# =============================================================================
echo "=============================================="
echo "Step 4: Generating metrics for comparison..."
echo "=============================================="

# Metrics for tumor-based reference
$CNVKIT metrics $OUTPUT_TUMOR_REF/*.cnr -s $OUTPUT_TUMOR_REF/*.cns \
    -o "$OUTPUT_TUMOR_REF/metrics_summary.tsv" 2>/dev/null || true

# Metrics for flat reference
$CNVKIT metrics $OUTPUT_FLAT_REF/*.cnr -s $OUTPUT_FLAT_REF/*.cns \
    -o "$OUTPUT_FLAT_REF/metrics_summary.tsv" 2>/dev/null || true

echo "=============================================="
echo "Pipeline completed!"
echo "=============================================="
echo ""
echo "Results:"
echo "  Tumor-based reference: $OUTPUT_TUMOR_REF/"
echo "  Flat reference:        $OUTPUT_FLAT_REF/"
echo ""
echo "Compare metrics:"
echo "  cat $OUTPUT_TUMOR_REF/metrics_summary.tsv"
echo "  cat $OUTPUT_FLAT_REF/metrics_summary.tsv"
echo ""
echo "Key metrics to compare:"
echo "  - spread: lower is better (less noise)"
echo "  - bivar:  lower is better (less variance)"
