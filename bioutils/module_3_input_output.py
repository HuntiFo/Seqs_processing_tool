import os


def fastq_transfom_dict(path: str) -> dict:
    counter = 4
    seqs = {}
    with open(path) as file:
        for line in file:
            edited_line = line.strip()
            if counter % 4 == 0:
                seqs[edited_line] = []
                new_line = edited_line
            elif counter % 4 == 1:
                seqs[new_line].append(edited_line)
            elif counter % 4 == 3:
                seqs[new_line].append(edited_line)
                seqs[new_line] = tuple(seqs[new_line])
            counter += 1
    return seqs


def dict_transform_fastq(
    filtered_reads: dict[str, tuple[str, str]], output_fastq: str
) -> None:
    if not os.path.isdir("filtered"):
        os.mkdir("filtered")
    if os.path.isfile("./filtered/" + output_fastq):
        return "[UPD]:file with this name already esists. Rename output file"
    with open("filtered" + os.path.sep + output_fastq, "w") as result_file:
        for seq_name, (sequence, quality) in filtered_reads.items():
            result_file.write(seq_name + "\n")
            result_file.write(sequence + "\n")
            result_file.write("+" + seq_name + "\n")
            result_file.write(quality + "\n")
