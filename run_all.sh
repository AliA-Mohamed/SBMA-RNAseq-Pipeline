#!/bin/bash
# Run all pipeline steps for both contrasts
set -e

CONFIG="config/config.yaml"
DATASET="SBMA_Fibroblast_Okada"
CONTRASTS=("purified_sbma_vs_control" "unpurified_sbma_vs_control")

for CONTRAST in "${CONTRASTS[@]}"; do
    echo ""
    echo "============================================================"
    echo "  Running pipeline for: $DATASET / $CONTRAST"
    echo "============================================================"

    for STEP in 02 03 04 05 06 07 08 09; do
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
