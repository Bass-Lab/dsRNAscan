#!/usr/bin/env python3
"""
Enhanced Output Module for dsRNA Overlap Analyzer
Provides user-friendly, colorful, and informative output formatting
"""

import sys
import os
from typing import Dict, Optional, List, Tuple
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Check if terminal supports colors
HAS_COLOR = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

# Color codes for terminal output
class Colors:
    HEADER = '\033[95m' if HAS_COLOR else ''
    BLUE = '\033[94m' if HAS_COLOR else ''
    CYAN = '\033[96m' if HAS_COLOR else ''
    GREEN = '\033[92m' if HAS_COLOR else ''
    WARNING = '\033[93m' if HAS_COLOR else ''
    FAIL = '\033[91m' if HAS_COLOR else ''
    ENDC = '\033[0m' if HAS_COLOR else ''
    BOLD = '\033[1m' if HAS_COLOR else ''
    UNDERLINE = '\033[4m' if HAS_COLOR else ''

def print_header(text: str, width: int = 80):
    """Print a formatted header"""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*width}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text.center(width)}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*width}{Colors.ENDC}")

def print_section(text: str, width: int = 80):
    """Print a section header"""
    print(f"\n{Colors.CYAN}{Colors.BOLD}{text}{Colors.ENDC}")
    print(f"{Colors.CYAN}{'-'*len(text)}{Colors.ENDC}")

def print_metric(label: str, value: str, color: str = None, indent: int = 2):
    """Print a formatted metric"""
    color_code = getattr(Colors, color.upper(), '') if color else ''
    print(f"{' '*indent}- {label}: {color_code}{Colors.BOLD}{value}{Colors.ENDC}")

def format_number(n: float, decimals: int = 1) -> str:
    """Format a number with thousands separator"""
    if isinstance(n, int) or n.is_integer():
        return f"{int(n):,}"
    return f"{n:,.{decimals}f}"

def get_enrichment_interpretation(fold_enrichment: float, p_value: float, z_score: float) -> Dict[str, str]:
    """Get biological interpretation of enrichment results"""
    interpretation = {
        'status': '',
        'strength': '',
        'confidence': '',
        'biological_meaning': '',
        'recommendation': ''
    }
    
    # Determine status
    if p_value >= 0.05:
        interpretation['status'] = 'NOT SIGNIFICANT'
        interpretation['color'] = 'WARNING'
    elif fold_enrichment > 1:
        interpretation['status'] = 'ENRICHED'
        interpretation['color'] = 'GREEN'
    else:
        interpretation['status'] = 'DEPLETED'
        interpretation['color'] = 'FAIL'
    
    # Determine strength
    if abs(z_score) < 2:
        interpretation['strength'] = 'weak'
    elif abs(z_score) < 5:
        interpretation['strength'] = 'moderate'
    elif abs(z_score) < 10:
        interpretation['strength'] = 'strong'
    else:
        interpretation['strength'] = 'very strong'
    
    # Determine confidence
    if p_value < 0.001:
        interpretation['confidence'] = 'extremely high confidence (p < 0.001)'
    elif p_value < 0.01:
        interpretation['confidence'] = 'high confidence (p < 0.01)'
    elif p_value < 0.05:
        interpretation['confidence'] = 'moderate confidence (p < 0.05)'
    else:
        interpretation['confidence'] = 'not confident (p ≥ 0.05)'
    
    # Biological meaning
    if p_value < 0.05:
        if fold_enrichment > 10:
            interpretation['biological_meaning'] = 'Strong functional association likely'
        elif fold_enrichment > 5:
            interpretation['biological_meaning'] = 'Moderate functional association suggested'
        elif fold_enrichment > 2:
            interpretation['biological_meaning'] = 'Weak functional association possible'
        elif fold_enrichment > 1:
            interpretation['biological_meaning'] = 'Marginal association detected'
        elif fold_enrichment < 0.5:
            interpretation['biological_meaning'] = 'Strong mutual exclusion likely'
        else:
            interpretation['biological_meaning'] = 'Moderate avoidance detected'
    else:
        interpretation['biological_meaning'] = 'No clear functional relationship'
    
    # Recommendations
    if p_value < 0.05 and fold_enrichment > 5:
        interpretation['recommendation'] = 'Consider detailed mechanistic studies'
    elif p_value < 0.05 and fold_enrichment > 2:
        interpretation['recommendation'] = 'Validate with additional datasets'
    elif p_value < 0.05:
        interpretation['recommendation'] = 'Interesting but needs more evidence'
    else:
        interpretation['recommendation'] = 'May need larger sample size or different approach'
    
    return interpretation

def print_enhanced_results(basic_results: Dict, enrichment_results: Optional[Dict] = None, 
                          args: Optional[object] = None):
    """Print enhanced, user-friendly results"""
    
    # Header
    print_header("dsRNA OVERLAP ANALYSIS RESULTS", 80)
    
    # Input information
    print_section("Input Information")
    print_metric("File", args.input_file if args else "Unknown")
    print_metric("Format", basic_results.get('format', 'Unknown'))
    print_metric("dsRNA subset", args.subset if args else "all")
    print_metric("Analysis type", "Strand-specific" if args and args.strand_specific else "Strand-agnostic")
    if args and args.min_overlap > 1:
        print_metric("Minimum overlap", f"{args.min_overlap} bp")
    
    # Dataset statistics
    print_section("Dataset Statistics")
    print_metric("Query features", format_number(basic_results['n_query_features']), 'BLUE')
    print_metric("dsRNA regions", format_number(basic_results['n_dsrna_regions']), 'BLUE')
    print_metric("dsRNA intervals", format_number(basic_results['n_dsrna_intervals']), 'BLUE')
    
    # Overlap results
    print_section("Overlap Results")
    
    query_pct = basic_results['query_overlap_pct']
    dsrna_pct = basic_results['dsrna_overlap_pct']
    
    # Color code based on percentage
    query_color = 'GREEN' if query_pct > 1 else 'WARNING' if query_pct > 0.1 else 'FAIL'
    dsrna_color = 'GREEN' if dsrna_pct > 1 else 'WARNING' if dsrna_pct > 0.1 else 'FAIL'
    
    print_metric("Query features overlapping dsRNA", 
                f"{format_number(basic_results['n_query_overlapping'])} ({query_pct:.1f}%)",
                query_color)
    print_metric("dsRNA intervals overlapping query", 
                f"{format_number(basic_results['n_dsrna_overlapping'])} ({dsrna_pct:.1f}%)",
                dsrna_color)
    print_metric("Reciprocal overlap score", 
                f"{basic_results['reciprocal_overlap_score']:.3f}")
    
    # Enrichment analysis if available
    if enrichment_results:
        print_section("Statistical Enrichment Analysis")
        
        # Sample sizes
        print_metric("Sample sizes", "")
        print_metric("Query features tested", format_number(basic_results['n_query_features']), indent=4)
        print_metric("dsRNA intervals tested", format_number(basic_results['n_dsrna_intervals']), indent=4)
        print_metric("Permutations performed", format_number(enrichment_results.get('n_permutations', 100)), indent=4)
        
        # Observed vs Expected
        print_metric("Overlap counts", "")
        print_metric("Observed overlaps", format_number(enrichment_results['observed_overlaps']), indent=4)
        print_metric("Expected overlaps", 
                    f"{enrichment_results['expected_overlaps']:.1f} ± {enrichment_results['control_std']:.1f}", indent=4)
        print_metric("Control range", 
                    f"{enrichment_results['control_min']} - {enrichment_results['control_max']}", indent=4)
        
        # Statistical metrics
        print_metric("Statistical tests", "")
        fold = enrichment_results['fold_enrichment']
        fold_color = 'GREEN' if fold > 2 else 'WARNING' if fold > 1 else 'FAIL'
        print_metric("Fold enrichment", f"{fold:.2f}x", fold_color, indent=4)
        
        z_score = enrichment_results['z_score']
        z_color = 'GREEN' if abs(z_score) > 3 else 'WARNING' if abs(z_score) > 2 else 'FAIL'
        print_metric("Z-score", f"{z_score:.2f}", z_color, indent=4)
        
        p_val = enrichment_results['p_value']
        p_color = 'GREEN' if p_val < 0.01 else 'WARNING' if p_val < 0.05 else 'FAIL'
        p_display = f"{p_val:.4f}" if p_val >= 0.0001 else f"{p_val:.2e}"
        print_metric("P-value (two-tailed)", p_display, p_color, indent=4)
        
        # Simple result statement
        if p_val < 0.05:
            if fold > 1:
                status_icon = "✅"
                status_text = f"SIGNIFICANT ENRICHMENT ({fold:.1f}x, p={p_display})"
                status_color = Colors.GREEN
            else:
                status_icon = "⚠️"
                status_text = f"SIGNIFICANT DEPLETION ({fold:.1f}x, p={p_display})"
                status_color = Colors.WARNING
        else:
            status_icon = "❌"
            status_text = f"NOT SIGNIFICANT (p={p_display})"
            status_color = Colors.FAIL
        
        print(f"\n  {status_icon} {status_color}{Colors.BOLD}{status_text}{Colors.ENDC}")
        
    # Summary
    print_section("Summary")
    
    if enrichment_results:
        # Quick summary table
        print("  Statistical summary:")
        print(f"    - N query features: {format_number(basic_results['n_query_features'])}")
        print(f"    - N overlapping: {format_number(basic_results['n_query_overlapping'])} ({basic_results['query_overlap_pct']:.1f}%)")
        if 'fold_enrichment' in enrichment_results:
            print(f"    - Fold change: {enrichment_results['fold_enrichment']:.2f}x")
            print(f"    - P-value: {enrichment_results['p_value']:.4f}")
            print(f"    - Significance level: {'***' if enrichment_results['p_value'] < 0.001 else '**' if enrichment_results['p_value'] < 0.01 else '*' if enrichment_results['p_value'] < 0.05 else 'ns'}")
    else:
        print(f"  - Total features analyzed: {format_number(basic_results['n_query_features'])}")
        print(f"  - Features overlapping dsRNA: {format_number(basic_results['n_query_overlapping'])} ({basic_results['query_overlap_pct']:.1f}%)")
        print(f"  - No statistical enrichment test performed (use --permutations N)")
    
    # Footer with timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n{Colors.CYAN}Analysis completed at {timestamp}{Colors.ENDC}")
    print(f"{Colors.CYAN}{'='*80}{Colors.ENDC}\n")

def create_enrichment_plot(enrichment_results: Dict, output_file: str = "enrichment_plot.png"):
    """Create a visual representation of enrichment analysis"""
    import matplotlib
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot 1: Observed vs Expected
    ax1 = axes[0]
    observed = enrichment_results['observed_overlaps']
    expected = enrichment_results['expected_overlaps']
    std = enrichment_results['control_std']
    
    bars = ax1.bar(['Expected', 'Observed'], [expected, observed], 
                   yerr=[std, 0], capsize=5,
                   color=['#808080', '#2ecc71' if observed > expected else '#e74c3c'],
                   edgecolor='black', linewidth=1)
    ax1.set_ylabel('Number of Overlaps', fontsize=11)
    ax1.set_title('Observed vs Expected Overlaps', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, max(observed, expected) * 1.2)
    
    # Add value labels on bars
    for bar, val in zip(bars, [expected, observed]):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.0f}', ha='center', va='bottom')
    
    # Plot 2: Effect Size and Significance
    ax2 = axes[1]
    fold = enrichment_results['fold_enrichment']
    p_val = enrichment_results['p_value']
    
    # Determine significance level
    if p_val < 0.001:
        sig_text = '***'
    elif p_val < 0.01:
        sig_text = '**'
    elif p_val < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'
    
    # Bar for fold enrichment
    bar_color = '#2ecc71' if fold > 1 and p_val < 0.05 else '#e74c3c' if fold < 1 and p_val < 0.05 else '#95a5a6'
    bars2 = ax2.bar(['Fold\nEnrichment'], [fold], color=bar_color, edgecolor='black', linewidth=1)
    
    # Add horizontal line at y=1
    ax2.axhline(y=1, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # Add value and significance
    ax2.text(0, fold + 0.1, f'{fold:.1f}x', ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax2.text(0, fold + 0.3, sig_text, ha='center', va='bottom', fontsize=14, fontweight='bold')
    
    ax2.set_ylabel('Fold Change', fontsize=11)
    ax2.set_title('Enrichment Effect Size', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, max(fold, 1) * 1.5)
    
    # Plot 3: Statistical Summary Table
    ax3 = axes[2]
    ax3.axis('off')
    
    z_score = enrichment_results['z_score']
    n_perms = enrichment_results.get('n_permutations', 100)
    
    # Table data
    table_data = [
        ['Metric', 'Value'],
        ['─' * 12, '─' * 12],
        ['Observed', f'{observed:,}'],
        ['Expected', f'{expected:.1f} ± {std:.1f}'],
        ['Fold change', f'{fold:.2f}x'],
        ['Z-score', f'{z_score:.2f}'],
        ['P-value', f'{p_val:.4f}'],
        ['Permutations', f'{n_perms:,}'],
        ['─' * 12, '─' * 12],
        ['Result', 'ENRICHED' if fold > 1 and p_val < 0.05 else 'DEPLETED' if fold < 1 and p_val < 0.05 else 'NOT SIG.']
    ]
    
    # Format as text
    text = ''
    for row in table_data:
        if row[0].startswith('─'):
            text += f'{row[0]:<12} {row[1]:<12}\n'
        else:
            text += f'{row[0]:<12} {row[1]:>12}\n'
    
    ax3.text(0.1, 0.5, text, fontsize=10, verticalalignment='center', 
             fontfamily='monospace', transform=ax3.transAxes)
    ax3.set_title('Statistical Summary', fontsize=12, fontweight='bold')
    
    plt.suptitle('dsRNA Overlap Enrichment Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\n{Colors.GREEN}📊 Enrichment plot saved to: {output_file}{Colors.ENDC}")
    
    return output_file

def print_subset_recommendation(n_features: int):
    """Recommend appropriate dsRNA subset based on input size"""
    print_section(" Subset Recommendation")
    
    if n_features < 100:
        print("  Your dataset is small. Recommended subsets:")
        print("    - 'conserved_high_conf' - Most stringent, best for small datasets")
        print("    - 'conserved' - Good balance of size and quality")
    elif n_features < 1000:
        print("  Your dataset is medium-sized. Recommended subsets:")
        print("    - 'conserved' - Good starting point")
        print("    - 'ml_structure' - ML-filtered, larger set")
        print("    - 'nonalu_high_conf' - Focus on non-repetitive elements")
    elif n_features < 10000:
        print("  Your dataset is large. Recommended subsets:")
        print("    - 'any_high_conf' - Any ML model high confidence")
        print("    - 'ml_gtex' - GTEx model predictions")
        print("    - Test both 'alu_high_conf' and 'nonalu_high_conf' separately")
    else:
        print("  Your dataset is very large. You can use:")
        print("    - 'all' - Complete dsRNA set (may be slow)")
        print("    - 'any_high_conf' - Good balance of size and quality")
        print("    - Consider strand-specific analysis for better resolution")