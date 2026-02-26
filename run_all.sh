#!/bin/bash
# Run all pipeline steps for configured contrasts
#
# Usage:
#   bash run_all.sh                    # Run steps 01-10 for all contrasts
#   bash run_all.sh --with-deseq2      # Include Step 00 (DESeq2 from raw counts)
#
# Step 00 (DESeq2) is run separately before the main loop because it requires
# its own arguments (--counts, --samples, --design, --contrast, --output).
# See scripts/00_deseq2_from_counts.R --help for details.

set -e

CONFIG="config/config.yaml"
DATASET="SBMA_Fibroblast_Okada"
CONTRASTS=("purified_sbma_vs_control" "unpurified_sbma_vs_control")

# Optional Step 00: DESeq2 from raw counts
if [[ "$1" == "--with-deseq2" ]]; then
    echo ""
    echo "============================================================"
    echo "  Step 00: DESeq2 from Raw Counts"
    echo "============================================================"
    echo "  Step 00 requires manual invocation with dataset-specific"
    echo "  arguments. Example:"
    echo ""
    echo "  Rscript scripts/00_deseq2_from_counts.R \\"
    echo "    --counts data/MyDataset/count_matrix.csv \\"
    echo "    --samples data/MyDataset/sample_info.csv \\"
    echo "    --design '~ condition' \\"
    echo "    --contrast 'condition,SBMA,Control' \\"
    echo "    --output data/MyDataset/deseq2_degs.csv"
    echo ""
    echo "  Then update config.yaml to point to the output file."
    echo "============================================================"
fi

for CONTRAST in "${CONTRASTS[@]}"; do
    echo ""
    echo "============================================================"
    echo "  Running pipeline for: $DATASET / $CONTRAST"
    echo "============================================================"

    for STEP in 01 02 03 04 05 06 07 08 09; do
        SCRIPT="scripts/${STEP}_*.R"
        SCRIPT_PATH=$(ls $SCRIPT 2>/dev/null | head -1)
        if [ -z "$SCRIPT_PATH" ]; then
            echo "WARNING: No script found for step $STEP"
            continue
        fi
        echo ""
        echo "--- Running step $STEP: $SCRIPT_PATH ---"
        Rscript "$SCRIPT_PATH" --config "$CONFIG" --dataset "$DATASET" --contrast "$CONTRAST" 2>&1 || {
            echo "WARNING: Step $STEP failed for $CONTRAST, continuing..."
        }
    done
done

# Step 10: Cross-comparison (runs across all contrasts)
echo ""
echo "============================================================"
echo "  Running Step 10: Cross-comparison"
echo "============================================================"
Rscript scripts/10_cross_comparison.R --config "$CONFIG" --dataset "$DATASET" --contrast purified_sbma_vs_control 2>&1 || {
    echo "WARNING: Step 10 failed, continuing..."
}

echo ""
echo "============================================================"
echo "  Pipeline complete!"
echo "============================================================"
