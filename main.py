from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from parser import parser
import sys
import logging

logging.basicConfig(
    filename="logs.log",
    filemode="w",
    format="{asctime} {levelname} +++ {message}",
    datefmt="%Y-%m-%d %H:%M:%S",
    style="{",
    level=logging.INFO,
    force=True,
    encoding="utf-8",
)


class BiologicalSequence:
    def __init__(self, data: str):
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        if isinstance(key, int):
            if key < 0 or key >= len(self.data):
                raise IndexError("Index out of range")
            return self.data[key]
        elif isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            return self.data[key]
        else:
            raise TypeError("Index must be an integer or slice")

    def __str__(self):
        return self.data

    def check_alphabet(self):
        amino_acids = {
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y",
        }
        data_upper = set(self.data.upper())
        if not data_upper.difference({"A", "T", "C", "G"}):
            return "DNA"
        elif not data_upper.difference({"A", "U", "C", "G"}):
            return "RNA"
        elif not data_upper.difference(amino_acids):
            return "Protein"
        else:
            raise ValueError("Uknown type of BiologicalSequence")


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, data=None):
        if self.__class__ is NucleicAcidSequence:
            raise NotImplementedError("This method is not implemented.")
        super().__init__(data)

    def reverse(self):
        return self.data[::-1]

    def complement(self):
        dic_complement = {
            "A": "T",
            "T": "A",
            "a": "t",
            "t": "a",
            "C": "G",
            "G": "C",
            "c": "g",
            "g": "c",
            "U": "A",
            "u": "a",
        }
        new_seq = "".join([dic_complement[n] for n in self.data])
        return new_seq

    def reverse_complement(self):
        return (self.complement())[::-1]


class DNASequence(NucleicAcidSequence):
    def __init__(self, data: str):
        super().__init__(data)
        self.data = data
        if super().check_alphabet() != "DNA":
            raise ValueError("It is not DNA")

    def transcribe(self):
        return (self.data).replace("T", "U").replace("t", "u")


class RNASequence(NucleicAcidSequence):
    def __init__(self, data: str):
        super().__init__(data)
        self.data = data
        if super().check_alphabet() != "RNA":
            raise ValueError("It is not RNA")


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, data: str):
        super().__init__(data)
        self.data = data
        if super().check_alphabet() != "Protein":
            raise ValueError("It is not Protein")

    def calculate_weight(self):
        amino_acid_weights = {
            "A": 89.09,  # Alanine
            "C": 121.16,  # Cysteine
            "D": 133.10,  # Aspartic Acid
            "E": 147.13,  # Glutamic Acid
            "F": 165.19,  # Phenylalanine
            "G": 75.07,  # Glycine
            "H": 155.16,  # Histidine
            "I": 131.17,  # Isoleucine
            "K": 146.19,  # Lysine
            "L": 131.17,  # Leucine
            "M": 149.21,  # Methionine
            "N": 132.12,  # Asparagine
            "P": 115.13,  # Proline
            "Q": 146.14,  # Glutamine
            "R": 174.20,  # Arginine
            "S": 105.09,  # Serine
            "T": 119.12,  # Threonine
            "V": 117.15,  # Valine
            "W": 204.23,  # Tryptophan
            "Y": 181.19,  # Tyrosine
        }
        weight = 0
        for amino_acid in self.data:
            weight += amino_acid_weights[amino_acid.upper()]
        return f"Your protein weight: {weight:.2f} g/mol"


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: float | tuple[float, float] = (0, 2 ** 32),
    quality_threshold: float = 0,
) -> None:
    if isinstance(gc_bounds, float | int):
        gc_bounds = [gc_bounds, 100]
    elif isinstance(gc_bounds, tuple):
        gc_bounds = [gc_bounds[0], gc_bounds[1]]

    if isinstance(length_bounds, float | int):
        length_bounds = [length_bounds, 2 ** 32]
    elif isinstance(length_bounds, tuple):
        length_bounds = [length_bounds[0], length_bounds[1]]

    with open(input_fastq, "r") as input_fastq, open(output_fastq, "w") as output_fastq:
        for record in SeqIO.parse(input_fastq, "fastq"):
            if min(record.letter_annotations["phred_quality"]) >= quality_threshold:
                if gc_bounds[0] <= gc_fraction(record.seq) * 100 <= gc_bounds[1]:
                    if length_bounds[0] <= len(record.seq) <= length_bounds[1]:
                        SeqIO.write(record, output_fastq, "fastq")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    input_fastq = args.input_fastq
    output_fastq = args.output_fastq
    # parse gc_bounds
    if args.gc_bounds is None:
        gc_bounds = (0, 100)
    else:
        if len(args.gc_bounds) == 1:
            gc_bounds = int(args.gc_bounds[0])
        elif len(args.gc_bounds) == 2:
            gc_bounds = tuple(filter(lambda x: float(x), args.gc_bounds))
            if gc_bounds[0] > gc_bounds[1]:
                gc_bounds = (gc_bounds[1], gc_bounds[0])
                logging.error(
                    f"Mistake the order of gc-bounds. New order is {gc_bounds}"
                )
        else:
            raise Exception

    # parse length_bounds
    if args.length_bounds is None:
        length_bounds = (0, 2 ** 32)
    else:
        if len(args.length_bounds) == 1:
            length_bounds = int(args.length_bounds[0])
        elif len(args.length_bounds) == 2:
            length_bounds = tuple(filter(lambda x: int(x), args.length_bounds))
        else:
            raise Exception
    # parse quality_threshold
    if args.quality_threshold is None:
        quality_threshold = 0
    else:
        quality_threshold = int(args.quality_threshold[0])

    filter_fastq(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold)

    logging.info(
        f"Fastq filtration with options:\n"
        f"  - input file: {input_fastq}\n"
        f"  - output file: {output_fastq}\n"
        f"  - gc_bounds: {gc_bounds}\n"
        f"  - length_bounds: {length_bounds}\n"
        f"  - quality_threshold: {quality_threshold}\n"
        "Performed successfully."
    )
