#!/bin/bash
#
# CNVkit Pipeline Script (HPC Optimized with GNU Parallel)
# Uses sample.txt to distinguish tumor vs normal samples
#

set -e

# =============================================================================
# Fix temp directory for container compatibility (Apptainer/Singularity)
# Force TMPDIR to current directory (container's /local is read-only)
# =============================================================================
export TMPDIR="$(pwd)/.tmp"
mkdir -p "$TMPDIR"

# =============================================================================
# Configuration
# =============================================================================
export FASTA="/insomnia001/depts/pmg/users/sc5006/reference/Homo_sapiens_assembly38.fasta"
export BAM_DIR="/insomnia001/depts/pmg/users/sc5006/data/bam_2019/final_bams"
SAMPLE_FILE="sample.txt"
export REF_DIR="ref"
export OUTPUT_DIR="cnvkit_output"

# Number of parallel jobs (you have 10 CPUs)
NJOBS=8

mkdir -p "$REF_DIR" "$OUTPUT_DIR"

# =============================================================================
# Parse sample.txt and create temporary files with sample lists
# =============================================================================
echo "Parsing sample file..."

# Create temp files for tumor and normal sample lists (in current directory for container compatibility)
TUMOR_LIST="$(pwd)/.tumor_samples.tmp"
NORMAL_LIST="$(pwd)/.normal_samples.tmp"
rm -f "$TUMOR_LIST" "$NORMAL_LIST"
trap "rm -f $TUMOR_LIST $NORMAL_LIST; rm -rf $TMPDIR" EXIT

# Parse sample.txt (handles tab, comma, or space delimiters)
while IFS=$'\t, ' read -r sample_id sample_type || [[ -n "$sample_id" ]]; do
    [[ -z "$sample_id" ]] && continue
    [[ "$sample_id" == "sample_id" ]] && continue
    [[ "$sample_id" =~ ^# ]] && continue

    sample_type_lower=$(echo "$sample_type" | tr '[:upper:]' '[:lower:]' | tr -d '\r')

    if [[ "$sample_type_lower" == "tumor" ]]; then
        echo "$sample_id" >> "$TUMOR_LIST"
    elif [[ "$sample_type_lower" == "normal" || "$sample_type_lower" == "blood" ]]; then
        echo "$sample_id" >> "$NORMAL_LIST"
    fi
done < "$SAMPLE_FILE"

N_TUMOR=$(wc -l < "$TUMOR_LIST")
N_NORMAL=$(wc -l < "$NORMAL_LIST")
echo "Found $N_TUMOR tumor samples and $N_NORMAL normal samples"

# =============================================================================
# Step 1-3: Generate access, target, and antitarget files (if needed)
# =============================================================================
[[ ! -f "$REF_DIR/access.hg38.bed" ]] && \
    cnvkit.py access "$FASTA" -o "$REF_DIR/access.hg38.bed"

[[ ! -f "$REF_DIR/Exome-Agilent_V6.target.bed" ]] && \
    cnvkit.py target "$REF_DIR/Exome-Agilent_V6.bed" \
        --annotate "$REF_DIR/refFlat.txt" --split \
        -o "$REF_DIR/Exome-Agilent_V6.target.bed"

[[ ! -f "$REF_DIR/Exome-Agilent_V6.antitarget.bed" ]] && \
    cnvkit.py antitarget "$REF_DIR/Exome-Agilent_V6.target.bed" \
        -g "$REF_DIR/access.hg38.bed" \
        -o "$REF_DIR/Exome-Agilent_V6.antitarget.bed"

# =============================================================================
# Step 4: Run coverage for ALL samples in parallel
# =============================================================================
echo "Running coverage analysis in parallel..."

run_coverage() {
    sample_id=$1
    bam_file="$BAM_DIR/${sample_id}_final.bam"

    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file" >&2
        return 1
    fi

    # Target coverage
    [[ ! -f "$OUTPUT_DIR/${sample_id}.targetcoverage.cnn" ]] && \
        cnvkit.py coverage "$bam_file" "$REF_DIR/Exome-Agilent_V6.target.bed" \
            -o "$OUTPUT_DIR/${sample_id}.targetcoverage.cnn"

    # Antitarget coverage
    [[ ! -f "$OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn" ]] && \
        cnvkit.py coverage "$bam_file" "$REF_DIR/Exome-Agilent_V6.antitarget.bed" \
            -o "$OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn"

    echo "Completed: $sample_id"
}
export -f run_coverage

# Combine all samples and run in parallel
cat "$TUMOR_LIST" "$NORMAL_LIST" | parallel -j "$NJOBS" run_coverage {}

# =============================================================================
# Step 5: Build reference from NORMAL samples only
# =============================================================================
echo "Building reference from normal samples..."

if [[ ! -f "$REF_DIR/reference_2019.cnn" ]]; then
    # Build the file list using normal samples
    NORMAL_CNN_ARGS=""
    while read -r sample_id; do
        NORMAL_CNN_ARGS="$NORMAL_CNN_ARGS $OUTPUT_DIR/${sample_id}.targetcoverage.cnn"
        NORMAL_CNN_ARGS="$NORMAL_CNN_ARGS $OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn"
    done < "$NORMAL_LIST"

    cnvkit.py reference $NORMAL_CNN_ARGS \
        --fasta "$FASTA" \
        -o "$REF_DIR/reference_2019.cnn"
fi

# =============================================================================
# Step 6: Fix and segment TUMOR samples only (in parallel)
# =============================================================================
echo "Running fix and segment on tumor samples..."

run_fix_segment() {
    sample_id=$1

    # Fix
    [[ ! -f "$OUTPUT_DIR/${sample_id}.cnr" ]] && \
        cnvkit.py fix \
            "$OUTPUT_DIR/${sample_id}.targetcoverage.cnn" \
            "$OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn" \
            "$REF_DIR/reference_2019.cnn" \
            -o "$OUTPUT_DIR/${sample_id}.cnr"

    # Segment
    [[ ! -f "$OUTPUT_DIR/${sample_id}.cns" ]] && \
        cnvkit.py segment \
            "$OUTPUT_DIR/${sample_id}.cnr" \
            -o "$OUTPUT_DIR/${sample_id}.cns"

    # Call
    [[ ! -f "$OUTPUT_DIR/${sample_id}.call.cns" ]] && \
        cnvkit.py call \
            "$OUTPUT_DIR/${sample_id}.cns" \
            -o "$OUTPUT_DIR/${sample_id}.call.cns"

    echo "Completed fix/segment: $sample_id"
}
export -f run_fix_segment

cat "$TUMOR_LIST" | parallel -j "$NJOBS" run_fix_segment {}

echo "Pipeline completed! Results in: $OUTPUT_DIR"
