import pandas as pd
import argparse
import json
from collections import defaultdict
try:
    from consensus_scaffold_config import AGP_FILES, AGP_NAMES, AGP_PREFIX, FAI_FILES
except ImportError:
    # Fallback if config file doesn't exist (for Snakemake usage)
    AGP_FILES = []
    AGP_NAMES = []
    AGP_PREFIX = []
    FAI_FILES = []

def parse_fai_file(fai_file):
    """Parse FAI file and return dictionary of contig -> length"""
    fai_dict = {}
    with open(fai_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            fai_dict[fields[0]] = int(fields[1])
    return fai_dict

def parse_agp_file(agp_file, target_prefix):
    """Parse AGP file and return dictionary of chromosome -> list of contigs (ordered)"""
    contig_dict = {}
    contig_order_dict = {}
    
    with open(agp_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[4] == "W":  # Only process W (whole genome) entries
                if fields[0].split('_')[-1] == "RagTag":
                    # take out the prefix
                    chromosome = fields[0].split('_')[0]
                    if target_prefix in ["YAO#1#GWHDOOG", "YAO#2#GWHDQZJ"]:
                        target_str = str(int(chromosome[len(target_prefix):]))
                        if target_str == "24":
                            chromosome = "chrM"
                        elif target_str == "23":
                            if target_prefix == "YAO#1#GWHDOOG":
                                chromosome = "chrY"
                            else:
                                chromosome = "chrX"
                        else:
                            chromosome = "chr" + target_str
                    else:
                        chromosome = chromosome[len(target_prefix):]
                    if contig_dict.get(chromosome):
                        contig_dict[chromosome].add(fields[5])
                        contig_order_dict[chromosome].append(fields[5])
                    else:
                        contig_dict[chromosome] = {fields[5]}
                        contig_order_dict[chromosome] = [fields[5]]
                else:
                    if contig_dict.get("unplaced"):
                        contig_dict["unplaced"].add(fields[5])
                    else:
                        contig_dict["unplaced"] = {fields[5]}
    
    return contig_dict, contig_order_dict

def create_consensus_contigs(all_lists, file_names, fai_file, list_chrom, min_occurrences=2, len_threshold=100000, dump_file='aggregated_contigs.csv'):
    """
    Create ordered consensus contigs with weighted scoring based on position information.
    
    Algorithm:
    1. For each contig, assign weights 0-100 based on its position in each reference
    2. Normalize weights by averaging across references that contain the contig
    3. Sort contigs by normalized weight to preserve order information
    4. Filter contigs that appear in >= min_occurrences samples OR have length >= len_threshold
    """
    if len(file_names) < 2:
        print("Need at least 2 samples for consensus building")
        return None
    
    #rest_samples = file_names[1:]  # Skip first sample, use rest for consensus
    rest_samples = file_names
    
    print(f"Creating consensus contigs from {rest_samples}")
    print("="*80)
    
    fai_dict = parse_fai_file(fai_file)
    
    # Get all chromosomes from all samples
    all_chromosomes = set()
    for name in file_names:
        all_chromosomes.update(all_lists[name].keys())
    
    f = open(dump_file, 'w')

    dict_recording_contigs = {}
    dict_contigs_voting = defaultdict(dict)
    
    for chrom in sorted(all_chromosomes):
        if chrom not in list_chrom:
            continue
        print(f"\nProcessing chromosome: {chrom}")
        
        # Collect all contigs for this chromosome across all samples
        all_contigs = set()
        for sample in rest_samples:
            if chrom in all_lists[sample]:
                all_contigs.update(all_lists[sample][chrom])
        
        # Count occurrences of each contig
        contig_counts = defaultdict(int)
        for sample in rest_samples:
            if chrom in all_lists[sample]:
                sample_contigs = all_lists[sample][chrom]
                for contig in sample_contigs:
                    contig_counts[contig] += 1
        
        # Filter contigs that appear in at least min_occurrences samples OR have length >= len_threshold
        consensus_contigs = {contig for contig, count in contig_counts.items() 
                           if count >= min_occurrences or fai_dict.get(contig, 0) >= len_threshold}
        for contig, count in contig_counts.items():
            dict_contigs_voting[contig][chrom] = count
        
        print(f"Consensus contigs (>= {min_occurrences} samples or >= {len_threshold} bp): {len(consensus_contigs)}")
        
        # Calculate weighted scores for each contig
        contig_weights = {}
        
        for contig in consensus_contigs:
            weights = []
            references_with_contig = []
            
            for sample in rest_samples:
                if chrom in all_lists[sample]:
                    sample_contigs = all_lists[sample][chrom]
                    if contig in sample_contigs:
                        # Find position of contig in this sample
                        try:
                            position = sample_contigs.index(contig)
                            # Calculate weight: 0-100 based on position
                            # Earlier positions get higher weights
                            total_contigs = len(sample_contigs)
                            if total_contigs > 1:
                                weight = 100 * (1 - position / (total_contigs - 1))
                            else:
                                weight = 100  # Single contig gets max weight
                            weights.append(weight)
                            references_with_contig.append(sample)
                        except ValueError:
                            # Contig not found (shouldn't happen given our filtering)
                            continue
            
            # Normalize weight by averaging across references that contain this contig
            if weights:
                normalized_weight = sum(weights) / len(weights)
                contig_weights[contig] = {
                    'normalized_weight': normalized_weight,
                    'occurrence_count': len(weights)
                }
        
        # Sort contigs by normalized weight (descending order)
        sorted_contigs = sorted(contig_weights.items(), 
                              key=lambda x: x[1]['normalized_weight'], 
                              reverse=True)
        
        # Write to CSV: chromosome,contig1,contig2,contig3,...
        dict_recording_contigs[chrom] = sorted_contigs
        """
        f.write(f"{chrom}")
        for contig, weight_info in sorted_contigs:
            f.write(f",{contig}")
        f.write("\n")
        """
        print(f"  Ordered {len(sorted_contigs)} contigs for {chrom}")
    
    dict_contig_assigned_chrom = defaultdict(set)
    for contig, chroms in dict_contigs_voting.items():
        if len(chroms) > 1:
            sorted_chroms = sorted(chroms.items(), key=lambda x: x[1], reverse=True)
            max_occurence = sorted_chroms[0][1]
            for chrom, occurence in sorted_chroms:
                if occurence > max_occurence/2:
                    dict_contig_assigned_chrom[contig].add(chrom)
                else:
                    break
        else:
            dict_contig_assigned_chrom[contig] = set([list(chroms.keys())[0]])

    for chrom, contigs in dict_recording_contigs.items():
        f.write(f"{chrom}")
        for contig, weight_info in contigs:
            if chrom in dict_contig_assigned_chrom[contig]:
                f.write(f",{contig}")
            else:
                pass
        f.write("\n")
    f.close()
    print(f"\nConsensus contigs saved to: {dump_file}")

def main():
    parser = argparse.ArgumentParser(description='Create consensus scaffold from multiple AGP files')
    parser.add_argument('--config', type=str, help='JSON config file with AGP files, names, prefixes, and FAI file')
    parser.add_argument('--agp-files', nargs='+', help='List of AGP file paths')
    parser.add_argument('--agp-names', nargs='+', help='List of names for each AGP file')
    parser.add_argument('--agp-prefixes', nargs='+', help='List of prefixes for each AGP file')
    parser.add_argument('--fai-file', type=str, help='FAI file path')
    parser.add_argument('-o', '--output', type=str, default='aggregated_contigs.csv', help='Output CSV file')
    parser.add_argument('--min-occurrences', type=int, default=2, help='Minimum occurrences for consensus')
    parser.add_argument('--len-threshold', type=int, default=100000, help='Length threshold for consensus')
    parser.add_argument('--autosome-only', action='store_true', help='Only use autosomes for consensus')
    parser.add_argument('--X', action='store_true', help='Include chrX for consensus')
    parser.add_argument('--Y', action='store_true', help='Include chrY for consensus')
    parser.add_argument('--M', action='store_true', help='Include chrM for consensus')
    
    args = parser.parse_args()
    flag_autosome = args.autosome_only
    
    # Load config from JSON file if provided, otherwise use command-line args
    if args.config:
        with open(args.config, 'r') as f:
            config = json.load(f)
        agp_files = config['agp_files']
        file_names = config['agp_names']
        agp_prefixes = config['agp_prefixes']
        fai_file = config['fai_file']
    elif args.agp_files and args.agp_names and args.agp_prefixes and args.fai_file:
        agp_files = args.agp_files
        file_names = args.agp_names
        agp_prefixes = args.agp_prefixes
        fai_file = args.fai_file
    else:
        # Fallback to config file (for backward compatibility)
        agp_files = AGP_FILES
        file_names = AGP_NAMES
        agp_prefixes = AGP_PREFIX
        fai_file = FAI_FILES[0] if FAI_FILES else None
        
        if not agp_files:
            parser.error("Either --config, or --agp-files/--agp-names/--agp-prefixes/--fai-file must be provided")
    
    if len(agp_files) != len(file_names) or len(agp_files) != len(agp_prefixes):
        parser.error("Number of AGP files, names, and prefixes must match")
    
    print(f"Processing {len(agp_files)} AGP files...")
    print("Files:", file_names)
    print("="*80)
    
    # Parse all AGP files
    all_dicts = {}
    all_lists = {}
    for agp_file, name, prefix in zip(agp_files, file_names, agp_prefixes):
        print(f"Parsing {name}...")
        all_dicts[name], all_lists[name] = parse_agp_file(agp_file, prefix)
    
    # Create consensus contigs
    list_chrom = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    if args.X:
        list_chrom.append("chrX")
    if args.Y:
        list_chrom.append("chrY")
    if args.M:
        list_chrom.append("chrM")
    create_consensus_contigs(
        all_lists, 
        file_names, 
        fai_file, 
        list_chrom=list_chrom,
        min_occurrences=args.min_occurrences, 
        len_threshold=args.len_threshold,
        dump_file=args.output,
    )
    
    print("\nDone!")

if __name__ == "__main__":
    main()

