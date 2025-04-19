import argparse

parser = argparse.ArgumentParser(
    prog="FASTQ FILTRATOR",
    description="This tool filtrate nucleotide reads from fastq file by gc_content, length and quality treshold.",
    epilog="Enjoy FASTQ FILTRATOR tool!",
)

parser.add_argument("input_fastq", type=str, help="Takes path to input fastq file")
parser.add_argument(
    "output_fastq", type=str, help="Takes path to output filtrated fastq file"
)
parser.add_argument(
    "-gc",
    "--gc_bounds",
    metavar="%%",
    type=float,
    nargs="+",
    help="Takes one or two floats as the limitation of gc content (in %%). Default from 0 to 100 percentages.",
)
parser.add_argument(
    "-lb",
    "--length_bounds",
    metavar="bp",
    type=int,
    nargs="+",
    help="Takes one or two intagers as the limitation of length of reads. Default from 0 to 2 ** 32 bp.",
)
parser.add_argument(
    "-qt",
    "--quality_threshold",
    type=int,
    nargs=1,
    help="Takes one intager as the limitation of read quality. Default is 0.",
)
