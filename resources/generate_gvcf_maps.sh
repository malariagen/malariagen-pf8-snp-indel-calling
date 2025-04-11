#!/bin/bash

# Define the pipeline directory
# PIPELINE_DIR="<INSERT PATH HERE>/pf9/"
GVCF_DIR="$PIPELINE_DIR/bams_gvcfs"
GVCF_MAP_DIR="$PIPELINE_DIR/gvcf_maps"

# Create the gvcf_map directory if it doesn't exist
mkdir -p "$GVCF_MAP_DIR"

# Define chromosome names
CHROMOSOMES=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "MIT" "API")

# Loop through each chromosome
for CHR in "${CHROMOSOMES[@]}"; do
    MAP_FILE="$GVCF_MAP_DIR/Pf3D7_${CHR}_v3.gvcf_map"
    echo "Generating $MAP_FILE"

    # Clear existing file or create a new one
    echo -n "" > "$MAP_FILE"

    # Find all sample directories
    for SAMPLE_DIR in "$GVCF_DIR"/*/output/gvcfs/*/; do
        SAMPLE_ID=$(basename "$SAMPLE_DIR")
        VCF_FILE=$(find "$SAMPLE_DIR" -type f -name "*Pf3D7_${CHR}_v3*.vcf.gz")
        
        if [[ -n "$VCF_FILE" ]]; then
            echo -e "$SAMPLE_ID\t$VCF_FILE" >> "$MAP_FILE"
        fi
    done

done

echo "GVCF map files generated successfully in $GVCF_MAP_DIR."
