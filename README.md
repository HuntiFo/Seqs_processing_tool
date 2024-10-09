# Seqs_processing_tool
This is a bioinformatic package for processing nucleotide sequences, which consists of two main functions `run_dna_rna_tools` and `filter_fastq`.

A) `run_dna_rna_tools` 
input a list of an arbitrary number of arguments with DNA or RNA sequences (str), as well as the name of the procedure to be performed (this is always the last argument, str). 
The list of actions:
1) `transcribe` — returns a transcribed sequnce
2) `reverse` — returns a reversed sequnce
3) `complement` — return a complementary sequence
4) `reverse_complement` — returns a reverse-complementary sequence
5) `gc_counter` — returns GC-composition in %
6) `nucacid_type` — returns the type of nucleic acid (RNA,DNA or RNA\DNA)
7) `is_aug_in_rna` — returns True/False, if AUG in RNA

After specified action on all the transmitted sequences, run_dna_rna_tools returns the result.

B) filter_fastq 
Input four arguments: seqs(dict), gc_bounds(float|tuple), length_bounds(float|tuple), quality_threshold(float).
1) `seqs` - dictionary, consisting of fastq reads.Key - str, the name of sequence. Value - a tuple of two lines: sequence of reads and quality of reads.
2) `gc_bounds` - str or tuple, the GC interval of the composition (in percent) for filtering (by default (0, 100)). If you pass a single number to the argument, it is assumed that this is the upper bound. Reads with GC composition higher than upper bound and lower than lower bound are discarded.
3) `length_bounds` - the length interval for filtering (by default it is (0, 2**32)). Reads with length longer than upper bound and shorter than lower bound are discarded.
4) `quality_threshold` - the threshold value of the average read quality for filtering (0 by       default) (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded.

Returns a dictionary with filtered reads by prescribed conditions. All the described intervals include both upper and lower bounds.
