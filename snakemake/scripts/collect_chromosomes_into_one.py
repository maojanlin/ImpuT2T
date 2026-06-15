import argparse



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input paths for the chromosomes to be aggregated")
    parser.add_argument("-n", "--names", nargs="+", required=True, help="Names for the chromosomes to be aggregated")
    parser.add_argument("-o", "--output", required=True, help="Output path for the aggregated genome FASTA")
    parser.add_argument(
        "-a",
        "--agp-input",
        nargs="+",
        required=False,
        help="Optional: input AGP paths corresponding to the input FASTA files (same order)",
    )
    parser.add_argument(
        "-g",
        "--agp-output",
        required=False,
        help="Optional: output path for the aggregated AGP",
    )
    args = parser.parse_args()

    input_paths = args.input
    names = args.names
    output_path = args.output

    # Aggregate FASTA sequences: prefix patch object names (column 1 in AGP).
    with open(output_path, "w") as f_o:
        for input_path, name in zip(input_paths, names):
            with open(input_path, "r") as f_i:
                for line in f_i:
                    if line.startswith(">"):
                        header = line.strip()[1:]
                        f_o.write(f">{name}_{header}\n")
                    else:
                        f_o.write(line.strip() + "\n")

    # Aggregate AGP: prefix object IDs (column 1) only; keep component IDs (column 6) as-is.
    if args.agp_input and args.agp_output:
        agp_input_paths = args.agp_input
        agp_output_path = args.agp_output

        with open(agp_output_path, "w") as f_agp_o:
            for agp_path, name in zip(agp_input_paths, names):
                with open(agp_path, "r") as f_agp_i:
                    for line in f_agp_i:
                        if line.startswith("#") or not line.strip():
                            continue

                        cols = line.rstrip("\n").split("\t")
                        if len(cols) < 6:
                            f_agp_o.write(line)
                            continue

                        cols[0] = f"{name}_{cols[0]}"
                        f_agp_o.write("\t".join(cols) + "\n")



if __name__ == "__main__":
    main()