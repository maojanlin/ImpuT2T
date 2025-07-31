#!/usr/bin/env python3
"""
Script to tabularize DNADiff report files into a TSV format.
Extracts key metrics from multiple .report files and outputs them as a table.
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Union


def parse_dnadiff_report(file_path: str) -> Dict[str, Union[int, str]]:
    """
    Parse a DNADiff report file and extract the specified metrics.
    
    Args:
        file_path: Path to the .report file
        
    Returns:
        Dictionary mapping metric names to their values
    """
    metrics = {}
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        return {}
    
    # Define patterns for each metric
    patterns = {
        'TotalSeqs': r'TotalSeqs\s+(\d+)\s+(\d+)',
        'AlignedBases(Ref)': r'AlignedBases\s+(\d+)\(([^)]+)\)\s+(\d+)\([^)]+\)',
        'AlignedBases(Qry)': r'AlignedBases\s+(\d+)\([^)]+\)\s+(\d+)\(([^)]+)\)',
        'Alignments_1-to-1(#Seg)': r'1-to-1\s+(\d+)\s+(\d+)',
        'Alignments_1-to-1(AvgID)': r'1-to-1\s+\d+\s+\d+\s*\n\s*TotalLength\s+\d+\s+\d+\s*\n\s*AvgLength\s+[0-9.]+\s+[0-9.]+\s*\n\s*AvgIdentity\s+([0-9.]+)\s+([0-9.]+)',
        'Alignments_M-to-M(#Seg)': r'M-to-M\s+(\d+)\s+(\d+)',
        'Alignments_M-to-M(AvgID)': r'M-to-M\s+\d+\s+\d+\s*\n\s*TotalLength\s+\d+\s+\d+\s*\n\s*AvgLength\s+[0-9.]+\s+[0-9.]+\s*\n\s*AvgIdentity\s+([0-9.]+)\s+([0-9.]+)',
        'Breakpoints(Ref)': r'Breakpoints\s+(\d+)\s+(\d+)',
        'Breakpoints(Qry)': r'Breakpoints\s+(\d+)\s+(\d+)',
        'Relocations(Ref)': r'Relocations\s+(\d+)\s+(\d+)',
        'Relocations(Qry)': r'Relocations\s+(\d+)\s+(\d+)',
        'Translocation(Ref)': r'Translocations\s+(\d+)\s+(\d+)',
        'Translocation(Qry)': r'Translocations\s+(\d+)\s+(\d+)',
        'Inversions(Ref)': r'Inversions\s+(\d+)\s+(\d+)',
        'Inversions(Qry)': r'Inversions\s+(\d+)\s+(\d+)',
        'Insertions(Ref)': r'Insertions\s+(\d+)\s+(\d+)',
        'Insertions(Qry)': r'Insertions\s+(\d+)\s+(\d+)',
        'TandemIns(Ref)': r'TandemIns\s+(\d+)\s+(\d+)',
        'TandemIns(Qry)': r'TandemIns\s+(\d+)\s+(\d+)',
        'TotalSNPs': r'TotalSNPs\s+(\d+)\s+(\d+)',
        'TotalGSNPs': r'TotalGSNPs\s+(\d+)\s+(\d+)',
        'TotalIndels': r'TotalIndels\s+(\d+)\s+(\d+)',
        'TotalGIndels': r'TotalGIndels\s+(\d+)\s+(\d+)'
    }
    
    for metric_name, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            # For most metrics, we want the second value (QRY column)
            # For TotalSeqs, we'll use the first value as it's the same for both
            if metric_name == 'TotalSeqs':
                metrics[metric_name] = int(match.group(1))
            elif metric_name == 'AlignedBases(Ref)':
                # Format: number (percentage)
                value = int(match.group(1))
                percentage = match.group(2)
                metrics[metric_name] = f"{value} ({percentage})"
            elif metric_name == 'AlignedBases(Qry)':
                # Format: number (percentage)
                value = int(match.group(2))
                percentage = match.group(3)
                metrics[metric_name] = f"{value} ({percentage})"
            elif metric_name in ['Alignments_1-to-1(#Seg)', 'Alignments_M-to-M(#Seg)']:
                # For segment counts, use the second value (QRY column)
                metrics[metric_name] = int(match.group(2))
            elif metric_name in ['Alignments_1-to-1(AvgID)', 'Alignments_M-to-M(AvgID)']:
                # For average identity, use the second value (QRY column) and format as percentage
                value = float(match.group(2))
                metrics[metric_name] = f"{value:.4f}%"
            elif metric_name in ['Breakpoints(Ref)', 'Relocations(Ref)', 'Translocation(Ref)', 'Inversions(Ref)', 'Insertions(Ref)', 'TandemIns(Ref)']:
                # For Ref metrics, use the first value (REF column)
                metrics[metric_name] = int(match.group(1))
            else:
                metrics[metric_name] = int(match.group(2))
        else:
            # If not found, set to 0
            metrics[metric_name] = 0
    
    return metrics


def output_tsv(metrics_data: Dict[str, Dict[str, Union[int, str]]], output_file: Optional[str] = None):
    """
    Output the metrics data in TSV format.
    
    Args:
        metrics_data: Dictionary mapping file names to their metrics
        output_file: Optional output file path. If None, prints to stdout
    """
    if not metrics_data:
        print("No data to output.", file=sys.stderr)
        return
    
    # Define the order of metrics
    metric_order = [
        'TotalSeqs',
        'AlignedBases(Ref)',
        'AlignedBases(Qry)',
        'Alignments_1-to-1(#Seg)',
        'Alignments_1-to-1(AvgID)',
        'Alignments_M-to-M(#Seg)',
        'Alignments_M-to-M(AvgID)',
        'Breakpoints(Ref)',
        'Breakpoints(Qry)',
        'Relocations(Ref)',
        'Relocations(Qry)',
        'Translocation(Ref)',
        'Translocation(Qry)',
        'Inversions(Ref)',
        'Inversions(Qry)',
        'Insertions(Ref)',
        'Insertions(Qry)',
        'TandemIns(Ref)',
        'TandemIns(Qry)',
        'TotalSNPs',
        'TotalGSNPs',
        'TotalIndels',
        'TotalGIndels'
    ]
    
    # Get all file names (headers)
    file_names = list(metrics_data.keys())
    
    # Prepare output
    output_lines = []
    
    # Header row
    header = ['Metric'] + file_names
    output_lines.append('\t'.join(header))
    
    # Data rows
    for metric in metric_order:
        row = [metric]
        for file_name in file_names:
            value = metrics_data[file_name].get(metric, 0)
            row.append(str(value))
        output_lines.append('\t'.join(row))
    
    # Output
    if output_file:
        try:
            with open(output_file, 'w') as f:
                f.write('\n'.join(output_lines) + '\n')
            print(f"Output written to {output_file}")
        except Exception as e:
            print(f"Error writing to {output_file}: {e}", file=sys.stderr)
            # Fall back to stdout
            print('\n'.join(output_lines))
    else:
        print('\n'.join(output_lines))


def main():
    parser = argparse.ArgumentParser(
        description='Parse DNADiff report files and output metrics in TSV format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s file1.report file2.report file3.report
  %(prog)s *.report -o output.tsv
  %(prog)s -i file1.report file2.report -o results.tsv
        """
    )
    
    parser.add_argument(
        'report_files',
        nargs='+',
        help='DNADiff report files to process (.report files)'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output file path (default: print to stdout)'
    )
    
    parser.add_argument(
        '--basename',
        action='store_true',
        help='Use basename of files as column headers instead of full paths'
    )
    
    args = parser.parse_args()
    
    # Validate input files
    valid_files = []
    for file_path in args.report_files:
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} does not exist, skipping.", file=sys.stderr)
            continue
        if not file_path.endswith('.report'):
            print(f"Warning: File {file_path} does not have .report extension.", file=sys.stderr)
        valid_files.append(file_path)
    
    if not valid_files:
        print("Error: No valid report files provided.", file=sys.stderr)
        sys.exit(1)
    
    # Parse all report files
    metrics_data = {}
    for file_path in valid_files:
        # Use basename or full path as key
        key = os.path.basename(file_path) if args.basename else file_path
        metrics = parse_dnadiff_report(file_path)
        if metrics:
            metrics_data[key] = metrics
        else:
            print(f"Warning: Could not parse {file_path}, skipping.", file=sys.stderr)
    
    if not metrics_data:
        print("Error: No valid data could be extracted from the provided files.", file=sys.stderr)
        sys.exit(1)
    
    # Output the results
    output_tsv(metrics_data, args.output)


if __name__ == '__main__':
    main()
