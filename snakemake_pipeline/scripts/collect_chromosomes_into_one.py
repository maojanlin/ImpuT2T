import argparse



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input paths for the chromosomes to be aggregated")
    parser.add_argument("-n", "--names", nargs="+", required=True, help="Names for the chromosomes to be aggregated")
    parser.add_argument("-o", "--output", required=True, help="Output path for the aggregated genome")
    args = parser.parse_args()

    input_paths = args.input
    names = args.names
    output_path = args.output

    f_o = open(output_path, "w")
    for input_path, name in zip(input_paths, names):
        with open(input_path, "r") as f_i:
            for line in f_i:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    f_o.write(f">{name}_{header}\n")
                else:
                    f_o.write(line.strip() + "\n")
        f_i.close()
    f_o.close()


if __name__ == "__main__":
    main()