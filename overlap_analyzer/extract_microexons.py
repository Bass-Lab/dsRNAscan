#!/usr/bin/env python3
"""
Extract microexons from GENCODE GFF3 annotation file.
Microexons are typically defined as exons ≤30bp, though some definitions use up to 51bp.
"""

import gzip
import sys
from collections import defaultdict
import argparse

def parse_gff_attributes(attr_string):
    """Parse GFF3 attribute string into a dictionary."""
    attrs = {}
    for attr in attr_string.strip().split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
    return attrs

def extract_microexons(gff_file, max_length=30, min_length=1):
    """
    Extract microexons from GENCODE GFF3 file.
    
    Args:
        gff_file: Path to GENCODE GFF3 file (can be gzipped)
        max_length: Maximum length for microexons (default 30bp)
        min_length: Minimum length for microexons (default 1bp)
    
    Returns:
        List of microexon records
    """
    microexons = []
    exon_count = 0
    
    # Open file (handle gzipped)
    if gff_file.endswith('.gz'):
        f = gzip.open(gff_file, 'rt')
    else:
        f = open(gff_file, 'r')
    
    print(f"Processing {gff_file}...", file=sys.stderr)
    print(f"Looking for exons between {min_length}-{max_length}bp", file=sys.stderr)
    
    for line in f:
        if line.startswith('#'):
            continue
            
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
            
        # Only process exon features
        if fields[2] != 'exon':
            continue
            
        exon_count += 1
        if exon_count % 100000 == 0:
            print(f"Processed {exon_count} exons...", file=sys.stderr)
        
        chrom = fields[0]
        source = fields[1]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes = parse_gff_attributes(fields[8])
        
        # Calculate exon length
        exon_length = end - start + 1
        
        # Check if it's a microexon
        if min_length <= exon_length <= max_length:
            # Extract relevant attributes
            gene_id = attributes.get('gene_id', 'NA')
            gene_name = attributes.get('gene_name', 'NA')
            gene_type = attributes.get('gene_type', 'NA')
            transcript_id = attributes.get('transcript_id', 'NA')
            transcript_type = attributes.get('transcript_type', 'NA')
            exon_number = attributes.get('exon_number', 'NA')
            exon_id = attributes.get('exon_id', 'NA')
            
            microexon = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'length': exon_length,
                'strand': strand,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'gene_type': gene_type,
                'transcript_id': transcript_id,
                'transcript_type': transcript_type,
                'exon_number': exon_number,
                'exon_id': exon_id,
                'source': source
            }
            microexons.append(microexon)
    
    f.close()
    
    print(f"\nTotal exons processed: {exon_count}", file=sys.stderr)
    print(f"Microexons found: {len(microexons)}", file=sys.stderr)
    
    return microexons

def write_bed(microexons, output_file):
    """Write microexons to BED format."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("#chrom\tstart\tend\tname\tscore\tstrand\tgene_name\tgene_type\tlength\n")
        
        for me in microexons:
            # BED format is 0-based, GFF is 1-based
            bed_start = me['start'] - 1
            bed_end = me['end']
            
            # Create a unique name for the microexon
            name = f"{me['gene_name']}_{me['exon_id']}_{me['length']}bp"
            
            # Use length as score for visualization
            score = me['length']
            
            # Write BED line with extended fields
            f.write(f"{me['chrom']}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{me['strand']}\t")
            f.write(f"{me['gene_name']}\t{me['gene_type']}\t{me['length']}\n")

def write_gff(microexons, output_file):
    """Write microexons to GFF3 format."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("##gff-version 3\n")
        f.write("#description: Microexons extracted from GENCODE annotation\n")
        
        for me in microexons:
            # Build attributes string
            attrs = []
            attrs.append(f"ID=microexon:{me['exon_id']}")
            attrs.append(f"gene_id={me['gene_id']}")
            attrs.append(f"gene_name={me['gene_name']}")
            attrs.append(f"gene_type={me['gene_type']}")
            attrs.append(f"transcript_id={me['transcript_id']}")
            attrs.append(f"transcript_type={me['transcript_type']}")
            attrs.append(f"exon_number={me['exon_number']}")
            attrs.append(f"exon_id={me['exon_id']}")
            attrs.append(f"microexon_length={me['length']}")
            
            attr_string = ';'.join(attrs)
            
            # Write GFF line
            f.write(f"{me['chrom']}\t{me['source']}\tmicroexon\t{me['start']}\t{me['end']}\t.\t{me['strand']}\t.\t{attr_string}\n")

def analyze_microexons(microexons):
    """Analyze and report statistics about microexons."""
    print("\n=== Microexon Statistics ===", file=sys.stderr)
    
    # Length distribution
    length_dist = defaultdict(int)
    for me in microexons:
        length_dist[me['length']] += 1
    
    print("\nLength distribution:", file=sys.stderr)
    for length in sorted(length_dist.keys()):
        print(f"  {length}bp: {length_dist[length]} exons", file=sys.stderr)
    
    # Gene type distribution
    gene_type_dist = defaultdict(int)
    for me in microexons:
        gene_type_dist[me['gene_type']] += 1
    
    print("\nTop gene types with microexons:", file=sys.stderr)
    for gene_type, count in sorted(gene_type_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"  {gene_type}: {count} microexons", file=sys.stderr)
    
    # Unique genes with microexons
    unique_genes = set(me['gene_name'] for me in microexons)
    print(f"\nUnique genes with microexons: {len(unique_genes)}", file=sys.stderr)
    
    # Strand distribution
    strand_dist = defaultdict(int)
    for me in microexons:
        strand_dist[me['strand']] += 1
    
    print("\nStrand distribution:", file=sys.stderr)
    for strand, count in strand_dist.items():
        print(f"  {strand}: {count} microexons", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Extract microexons from GENCODE GFF3 annotation')
    parser.add_argument('--gff', '-g', required=True,
                        help='Path to GENCODE GFF3 file (can be gzipped)')
    parser.add_argument('--max-length', '-m', type=int, default=30,
                        help='Maximum length for microexons (default: 30bp)')
    parser.add_argument('--min-length', '-n', type=int, default=1,
                        help='Minimum length for microexons (default: 1bp)')
    parser.add_argument('--output-bed', '-b',
                        help='Output BED file path')
    parser.add_argument('--output-gff', '-o',
                        help='Output GFF3 file path')
    parser.add_argument('--stats', '-s', action='store_true',
                        help='Print statistics about microexons')
    
    args = parser.parse_args()
    
    # Extract microexons
    microexons = extract_microexons(args.gff, args.max_length, args.min_length)
    
    # Write output files
    if args.output_bed:
        write_bed(microexons, args.output_bed)
        print(f"Wrote BED file: {args.output_bed}", file=sys.stderr)
    
    if args.output_gff:
        write_gff(microexons, args.output_gff)
        print(f"Wrote GFF file: {args.output_gff}", file=sys.stderr)
    
    # Print statistics if requested
    if args.stats:
        analyze_microexons(microexons)
    
    # If no output file specified, write BED to stdout
    if not args.output_bed and not args.output_gff:
        for me in microexons:
            bed_start = me['start'] - 1
            name = f"{me['gene_name']}_{me['exon_id']}_{me['length']}bp"
            print(f"{me['chrom']}\t{bed_start}\t{me['end']}\t{name}\t{me['length']}\t{me['strand']}")

if __name__ == '__main__':
    main()