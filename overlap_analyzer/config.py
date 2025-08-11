"""
Configuration for dsRNA overlap analyzer within dsRNAscan
"""
import os
from pathlib import Path

def find_dsrna_data():
    """Find dsRNA data file in various locations"""
    
    # Possible locations relative to overlap_analyzer directory
    possible_paths = [
        # Minimal parquet in data directory (primary)
        Path(__file__).parent / "data" / "dsrna_predictions_minimal.parquet",
        
        # Full parquet if available
        Path(__file__).parent / "data" / "dsrna_predictions.parquet",
        
        # In parent dsRNAscan output
        Path(__file__).parent.parent / "output" / "dsrna_predictions.parquet",
        
        # User specified in environment variable
        os.environ.get("DSRNA_DATA_PATH", ""),
        
        # Common output names from dsRNAscan
        Path(__file__).parent.parent / "dsRNAscan_output.parquet",
        Path(__file__).parent.parent / "predictions.parquet",
    ]
    
    for path in possible_paths:
        if path and Path(path).exists():
            return str(path)
    
    return None

# Default data path
DEFAULT_DSRNA_FILE = find_dsrna_data()
