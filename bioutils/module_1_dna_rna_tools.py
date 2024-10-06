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


def transcribe(seq: str) -> str:
    seq = seq.replace("T", "U")
    seq = seq.replace("t", "u")
    return seq


def reverse(seq: str) -> str:
    return seq[::-1]


def complement(seq: str) -> str:
    new_seq = "".join([dic_complement[n] for n in seq])
    return new_seq


def reverse_complement(seq: str) -> str:
    return reverse(complement(seq))


def gc_counter(seq: str) -> float:
    gc_count = (seq.lower()).count("g") + (seq.lower()).count("c")
    gc_percents_seq = (gc_count / len(seq)) * 100
    return gc_percents_seq


def nucacid_type(seq: str) -> str:
    if "t" in seq.lower():
        return "DNA"
    elif "u" in seq.lower():
        return "RNA"
    else:
        return "RNA/DNA"


def is_aug_in_rna(seq: str) -> bool:
    if nucacid_type(seq) == "RNA":
        return seq.lower().count("aug") != 0
    return False
