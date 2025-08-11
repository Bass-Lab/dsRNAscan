#!/bin/bash
# Download dsRNA predictions data for overlap analyzer

echo "=========================================="
echo "dsRNA Data Download Helper"
echo "=========================================="
echo ""
echo "The full dsRNA predictions dataset is available from:"
echo ""
echo "1. Paper Supplementary Data:"
echo "   https://doi.org/10.1101/2025.01.24.634786"
echo ""
echo "2. Zenodo (when available):"
echo "   [DOI will be added upon publication]"
echo ""
echo "3. Direct download (example with wget):"
echo "   # Replace URL with actual data location"
echo "   wget -O data/dsrna_predictions.parquet [URL]"
echo ""
echo "After downloading, prepare subset files with:"
echo "   python prepare_data.py --input data/dsrna_predictions.parquet"
echo ""
echo "=========================================="

# Check if we can find the file locally
LOCAL_FILE="/Users/ryanandrews/Bioinformatics/20250619.df_with_normalized_predictions.parquet"

if [ -f "$LOCAL_FILE" ]; then
    echo ""
    echo "Found local dsRNA file at:"
    echo "  $LOCAL_FILE"
    echo ""
    read -p "Would you like to use this file? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Preparing data subsets..."
        python prepare_data.py --input "$LOCAL_FILE"
    fi
else
    echo ""
    echo "To use a local file, run:"
    echo "  python prepare_data.py --input /path/to/your/dsrna_predictions.parquet"
fi