def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> None:
    """
    Description of convert_multiline_fasta_to_oneline

    :param input_fasta: path to fasta file.
    :param output_fast: name of the new file.
    :return: None
    """
    with open(input_fasta) as fasta_file:
        with open(output_fasta, "w") as multiline_fasta_file:
            multiline_fasta_file.write(fasta_file.readline())
            for line in fasta_file:
                if line.startswith(">"):
                    multiline_fasta_file.write("\n" + line)
                else:
                    multiline_fasta_file.write(line.strip())


def parse_blast_output(input_fasta: str, output_fast: str) -> None:
    """
    Description of function parse_blast_output

    :param input_fasta: path to blast_result file.
    :param output_fast: name of the new file.
    :return: None
    """
    with open(input_fasta) as input_file:
        with open(output_fast, "w") as output_file:
            descriptions = []
            for line in input_file:
                if "Description" in line:
                    descriptions.append(
                        (input_file.readline().split("    ")[0]).split("...")[0]
                    )
                    descriptions.sort()
            for name in descriptions:
                output_file.write(name + "\n")
