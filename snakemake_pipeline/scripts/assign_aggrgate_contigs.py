import argparse
from pathlib import Path

def parse_fasta_file(fasta_file):
    dict_fasta = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_id = line.strip().split('>')[1]
                dict_fasta[contig_id] = ''
            else:
                dict_fasta[contig_id] += line.strip()
    return dict_fasta

def write_contigs_to_fasta(contigs_df, output_path, dict_fasta):
    """
    with open(output_path, 'w') as f:
        for contig_id in contigs_df['contig_id']:
            f.write(f">{contig_id}\n")
            f.write(f"{dict_fasta[contig_id]}\n")
    """
    for chromosome, contigs in contigs_df.items():
        f = open(f"{output_path}.{chromosome}.assigned.fa", 'w')
        for idx, contig_id in enumerate(contigs):
            f.write(f">{contig_id}_{chromosome}_impuT2T_{str(idx)}\n")
            f.write(f"{dict_fasta[contig_id]}\n")
        f.close()

def parse_contig_csv(contigs_csv):
    dict_chromosomes = {}
    with open(contigs_csv, 'r') as f:
        for line in f:
            fields = line.strip().split(',')
            chromosome = fields[0]
            contigs = fields[1:]
            dict_chromosomes[chromosome] = contigs
    return dict_chromosomes

def main():
    parser = argparse.ArgumentParser(description='Assign aggregate contigs to chromosomes')
    parser.add_argument('--contigs_file', type=str, required=True, help='Path to the contigs file')
    parser.add_argument('--contigs_csv', type=str, required=True, help='Path to the contigs dataframe')
    parser.add_argument('--output_path', type=str, default='./contigs_assignment/target', help='Path to the output file [./contigs_assignment/target]')
    args = parser.parse_args()

    dict_fasta = parse_fasta_file(args.contigs_file)
    dict_chromosomes = parse_contig_csv(args.contigs_csv)
    output_path = Path(args.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    write_contigs_to_fasta(dict_chromosomes, output_path, dict_fasta)


if __name__ == "__main__":
    main()
    



