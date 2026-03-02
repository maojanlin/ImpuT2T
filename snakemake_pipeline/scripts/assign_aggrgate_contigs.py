import argparse
from pathlib import Path

def parse_fasta_file(fasta_file):
    dict_fasta = {}
    contig_id = None
    contig_sequence = ''
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if contig_id is not None:
                    dict_fasta[contig_id] = contig_sequence
                contig_id = line.strip().split('>')[1]
                contig_sequence = ''
            else:
                contig_sequence += line.strip()
        dict_fasta[contig_id] = contig_sequence
    return dict_fasta

def write_contigs_to_fasta(contigs_df, output_path, dict_fasta):
    """
    with open(output_path, 'w') as f:
        for contig_id in contigs_df['contig_id']:
            f.write(f">{contig_id}\n")
            f.write(f"{dict_fasta[contig_id]}\n")
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for chromosome, contigs in contigs_df.items():
        out_file = output_path / f"{chromosome}.assigned.fa"
        with open(out_file, 'w') as f:
            for idx, contig_id in enumerate(contigs):
                f.write(f">{contig_id}_{chromosome}_impuT2T_{str(idx)}\n")
                f.write(f"{dict_fasta[contig_id]}\n")

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

    write_contigs_to_fasta(dict_chromosomes, args.output_path, dict_fasta)


if __name__ == "__main__":
    main()
    



