#!/bin/bash
#
# CNVkit Pipeline Script
# Uses sample.txt to distinguish tumor vs normal samples
#

set -e  # Exit on error

# =============================================================================
# Configuration
# =============================================================================
export FASTA="/insomnia001/depts/pmg/users/sc5006/reference/Homo_sapiens_assembly38.fasta"
BAM_DIR="/insomnia001/depts/pmg/users/sc5006/data/bam_2019/final_bams"
SAMPLE_FILE="sample.txt"  # Two columns: sample_id, Tumor/Normal
REF_DIR="ref"
OUTPUT_DIR="cnvkit_output"

# Number of parallel processes (adjust based on your resources)
NPROCS=10

# =============================================================================
# Setup directories
# =============================================================================
mkdir -p "$REF_DIR"
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# Step 1: Generate access file (if not already done)
# =============================================================================
if [[ ! -f "$REF_DIR/access.hg38.bed" ]]; then
    echo "Generating access file..."
    cnvkit.py access "$FASTA" -o "$REF_DIR/access.hg38.bed"
fi

# =============================================================================
# Step 2: Autobin (optional - to determine optimal bin sizes)
# You can skip this if you already know your bin sizes
# =============================================================================
# cnvkit.py autobin "$BAM_DIR"/*.bam \
#     -t "$REF_DIR/Exome-Agilent_V6.bed" \
#     -g "$REF_DIR/access.hg38.bed"

# =============================================================================
# Step 3: Generate target and antitarget BED files
# =============================================================================
if [[ ! -f "$REF_DIR/Exome-Agilent_V6.target.bed" ]]; then
    echo "Generating target BED..."
    cnvkit.py target "$REF_DIR/Exome-Agilent_V6.bed" \
        --annotate "$REF_DIR/refFlat.txt" \
        --split \
        -o "$REF_DIR/Exome-Agilent_V6.target.bed"
fi

if [[ ! -f "$REF_DIR/Exome-Agilent_V6.antitarget.bed" ]]; then
    echo "Generating antitarget BED..."
    cnvkit.py antitarget "$REF_DIR/Exome-Agilent_V6.target.bed" \
        -g "$REF_DIR/access.hg38.bed" \
        -o "$REF_DIR/Exome-Agilent_V6.antitarget.bed"
fi

# =============================================================================
# Step 4: Parse sample.txt and create tumor/normal lists
# =============================================================================
echo "Parsing sample file..."

# Arrays to store sample IDs
declare -a TUMOR_SAMPLES
declare -a NORMAL_SAMPLES

# Read sample.txt (skip header if present)
while IFS=$'\t, ' read -r sample_id sample_type || [[ -n "$sample_id" ]]; do
    # Skip empty lines and header
    [[ -z "$sample_id" ]] && continue
    [[ "$sample_id" == "sample_id" ]] && continue
    [[ "$sample_id" =~ ^# ]] && continue

    # Convert to lowercase for comparison
    sample_type_lower=$(echo "$sample_type" | tr '[:upper:]' '[:lower:]')

    if [[ "$sample_type_lower" == "tumor" ]]; then
        TUMOR_SAMPLES+=("$sample_id")
    elif [[ "$sample_type_lower" == "normal" || "$sample_type_lower" == "blood" ]]; then
        NORMAL_SAMPLES+=("$sample_id")
    else
        echo "Warning: Unknown sample type '$sample_type' for sample '$sample_id'"
    fi
done < "$SAMPLE_FILE"

echo "Found ${#TUMOR_SAMPLES[@]} tumor samples and ${#NORMAL_SAMPLES[@]} normal samples"

# =============================================================================
# Step 5: Run coverage for ALL samples (both tumor and normal)
# =============================================================================
echo "Running coverage analysis..."

# Function to run coverage for a single sample
run_coverage() {
    local sample_id=$1
    local bam_file="$BAM_DIR/${sample_id}_final.bam"

    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file"
        return 1
    fi

    # Target coverage
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.targetcoverage.cnn" ]]; then
        echo "  Processing target coverage: $sample_id"
        cnvkit.py coverage "$bam_file" \
            "$REF_DIR/Exome-Agilent_V6.target.bed" \
            -o "$OUTPUT_DIR/${sample_id}.targetcoverage.cnn"
    fi

    # Antitarget coverage
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn" ]]; then
        echo "  Processing antitarget coverage: $sample_id"
        cnvkit.py coverage "$bam_file" \
            "$REF_DIR/Exome-Agilent_V6.antitarget.bed" \
            -o "$OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn"
    fi
}

export -f run_coverage
export BAM_DIR REF_DIR OUTPUT_DIR

# Process all samples (tumor + normal)
ALL_SAMPLES=("${TUMOR_SAMPLES[@]}" "${NORMAL_SAMPLES[@]}")

# Option 1: Sequential processing
for sample_id in "${ALL_SAMPLES[@]}"; do
    run_coverage "$sample_id"
done

# Option 2: Parallel processing (uncomment to use)
# printf '%s\n' "${ALL_SAMPLES[@]}" | xargs -P "$NPROCS" -I {} bash -c 'run_coverage "$@"' _ {}

# =============================================================================
# Step 6: Build reference from NORMAL samples only
# =============================================================================
echo "Building reference from normal samples..."

if [[ ! -f "$REF_DIR/reference_2019.cnn" ]]; then
    # Build list of normal coverage files
    NORMAL_CNN_FILES=""
    for sample_id in "${NORMAL_SAMPLES[@]}"; do
        NORMAL_CNN_FILES="$NORMAL_CNN_FILES $OUTPUT_DIR/${sample_id}.targetcoverage.cnn"
        NORMAL_CNN_FILES="$NORMAL_CNN_FILES $OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn"
    done

    cnvkit.py reference $NORMAL_CNN_FILES \
        --fasta "$FASTA" \
        -o "$REF_DIR/reference_2019.cnn"
fi

# =============================================================================
# Step 7: Fix and segment TUMOR samples only
# =============================================================================
echo "Running fix and segment on tumor samples..."

for sample_id in "${TUMOR_SAMPLES[@]}"; do
    echo "  Processing tumor sample: $sample_id"

    # Fix
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.cnr" ]]; then
        cnvkit.py fix \
            "$OUTPUT_DIR/${sample_id}.targetcoverage.cnn" \
            "$OUTPUT_DIR/${sample_id}.antitargetcoverage.cnn" \
            "$REF_DIR/reference_2019.cnn" \
            -o "$OUTPUT_DIR/${sample_id}.cnr"
    fi

    # Segment
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.cns" ]]; then
        cnvkit.py segment \
            "$OUTPUT_DIR/${sample_id}.cnr" \
            -o "$OUTPUT_DIR/${sample_id}.cns"
    fi
done

# =============================================================================
# Step 8: Optional - Generate plots and calls for tumor samples
# =============================================================================
echo "Generating scatter plots and calls..."

for sample_id in "${TUMOR_SAMPLES[@]}"; do
    # Scatter plot
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.scatter.png" ]]; then
        cnvkit.py scatter \
            "$OUTPUT_DIR/${sample_id}.cnr" \
            -s "$OUTPUT_DIR/${sample_id}.cns" \
            -o "$OUTPUT_DIR/${sample_id}.scatter.png"
    fi

    # Diagram plot
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.diagram.pdf" ]]; then
        cnvkit.py diagram \
            "$OUTPUT_DIR/${sample_id}.cnr" \
            -s "$OUTPUT_DIR/${sample_id}.cns" \
            -o "$OUTPUT_DIR/${sample_id}.diagram.pdf"
    fi

    # Call copy number
    if [[ ! -f "$OUTPUT_DIR/${sample_id}.call.cns" ]]; then
        cnvkit.py call \
            "$OUTPUT_DIR/${sample_id}.cns" \
            -o "$OUTPUT_DIR/${sample_id}.call.cns"
    fi
done

echo "CNVkit pipeline completed!"
echo "Results are in: $OUTPUT_DIR"
