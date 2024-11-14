from bioutils.module_1_dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
    gc_counter,
    nucacid_type,
    is_aug_in_rna,
)

from bioutils.module_2_fastq import is_gc_fit, is_length_fit, is_quality_fit
from bioutils.module_3_input_output import fastq_transfom_dict, dict_transform_fastq


def run_dna_rna_tools(*args: str) -> str | list:
    """
    Description of function run_dna_rna_tools

    :param args: list of DNA/RNA sequences. Last element of this list is action
    :return: processed seqs made by action (if one seq -> reurn str,
    else -> return list of str)
    """
    *seqs, action = args
    result = []
    for seq in seqs:
        if "u" in seq.lower() and "t" in seq.lower():
            return "[UPD] The sequence contains both U and T at the same time"
        if action == "transcribe":
            result.append(transcribe(seq))
        elif action == "reverse":
            result.append(reverse(seq))
        elif action == "complement":
            result.append(complement(seq))
        elif action == "reverse_complement":
            result.append(reverse_complement(seq))
        elif action == "gc_counter":
            result.append(gc_counter(seq))
        elif action == "nucacid_type":
            result.append(nucacid_type(seq))
        elif action == "is_aug_in_rna":
            result.append(is_aug_in_rna(seq))
    if len(result) == 1:
        return result[0]
    return result


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: float | tuple[float, float] = (0, 2 ** 32),
    quality_threshold: float = 0,
) -> dict:
    """
    Description of function run_dna_rna_tools

    :param input_fastq: path to fastq file.
    :param output_fastq: new name of filtrated fastq file.
    :param gc_bounds: float|tuple
        The GC interval of the composition (in percent)
        for filtering (by default is (0, 100))
            float - define the upper bound
            tuple - define the interval
    :param length_bounds: float|tuple
        The length interval for filtering (by default it is (0, 2**32))
            float - define the upper bound
            tuple - define the interval
    :param quality_threshold: float
        The threshold value of the average read quality for filtering
        (by default 0 (phred33 scale))

    :return: processed seqs made by action (if one seq -> reurn str,
    else -> return list of str)
    """
    seqs = fastq_transfom_dict(input_fastq)
    filtered_reads = {}
    for seq_name, seq_data in seqs.items():
        seq, qulity = seq_data
        if (
            is_gc_fit(seq, gc_bounds)
            and is_length_fit(seq, length_bounds)
            and is_quality_fit(qulity, quality_threshold)
        ):
            filtered_reads[seq_name] = seq_data
    return dict_transform_fastq(filtered_reads, output_fastq)
