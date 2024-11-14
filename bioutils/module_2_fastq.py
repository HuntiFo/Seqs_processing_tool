from bioutils.module_1_dna_rna_tools import gc_counter


def is_gc_fit(seq: str, gc_bounds: float | tuple[float, float]) -> bool:
    percent_gc = gc_counter(seq)
    if isinstance(gc_bounds, tuple):
        return gc_bounds[0] <= percent_gc <= gc_bounds[1]
    return percent_gc <= gc_bounds


def is_length_fit(
    seq: str,
    length_bounds: float | tuple[float, float],
) -> bool:
    seq_lengh = len(seq)
    if isinstance(length_bounds, tuple):
        return length_bounds[0] <= seq_lengh <= length_bounds[1]
    return seq_lengh <= length_bounds


def is_quality_fit(qulity_seq: str, quality_threshold: float) -> bool:
    qulity_sum = 0
    for character in qulity_seq:
        qulity_sum += ord(character) - 33
    return qulity_sum / len(qulity_seq) >= quality_threshold
