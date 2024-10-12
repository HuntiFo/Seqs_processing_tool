# Seqs_processing_tool

This is a bioinformatic package for processing nucleotide sequences, which consists of two main functions `run_dna_rna_tools` and `filter_fastq`.

## Installation

To install Seqs_processing_tool you need to download the repository from GitHub

```
git clone https://github.com/HuntiFo/Seqs_processing_tool.git
cd Seqs_processing_tool
```

New directory Seqs_processing_tool contains:
1) the main script `Seqs_processing_tool.py` (with two fucntions `run_dna_rna_tools` and `filter_fastq`)
2) bio_files_processor.py 
README.md and directory bioutils with three modules: module_1_dna_rna_tools, module_2_fastq, module_3_input_ouput.

## Running instructions

### Main function

To run run_dna_rna_tools open the main file and write list of sequences and call the fuction run_dna_rna_tools(). 
A) `run_dna_rna_tools` 
input a list of an arbitrary number of arguments with DNA or RNA sequences (str), as well as the name of the procedure to be performed (this is always the last argument, str). 
As the last element of the list write the next actions:
1) `transcribe` — returns a transcribed sequnce
2) `reverse` — returns a reversed sequnce
3) `complement` — return a complementary sequence
4) `reverse_complement` — returns a reverse-complementary sequence
5) `gc_counter` — returns GC-composition in %
6) `nucacid_type` — returns the type of nucleic acid (RNA,DNA or RNA\DNA)
7) `is_aug_in_rna` — returns True/False, if AUG in RNA
After specified action on all the transmitted sequences, run_dna_rna_tools returns the result.

**Example:**
```
run_dna_rna_tools(['read1','read2',`transcribe`])
```

B) `filter_fastq`
To run `filter_fastq` open the main file, write the path to the fastq file, the name to the new file and next additional options. Call the function. New file will be created in the new directory filtered.
The function have five arguments: input_fastq(str),output_fastq(str), gc_bounds(float|tuple), length_bounds(float|tuple), quality_threshold(float).
1) `input_fastq` - str, the path to fastq reads.
2) `output_fastq` - str, the name of a new file with selected fastq reads.
3) `gc_bounds` - str or tuple, the GC interval of the composition (in percent) for filtering (by default (0, 100)). If you pass a single number to the argument, it is assumed that this is the upper bound. Reads with GC composition higher than upper bound and lower than lower bound are discarded.
3) `length_bounds` - the length interval for filtering (by default it is (0, 2**32)). Reads with length longer than upper bound and shorter than lower bound are discarded.
4) `quality_threshold` - the threshold value of the average read quality for filtering (0 by       default) (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded.
All the described intervals include both upper and lower bounds.

 **Example:**
```
filter_fastq(path_to_fastq, name_of_ouput_file, gc_bounds = (10,90), length_bounds = (30,90), quality_threshold = 5)
```
### Additional functions
There are two functions in bio_files_processor.py: convert_multiline_fasta_to_oneline and parse_blast_output. 
The `convert_multiline_fasta_to_oneline` convert multiple reads from fasta file in one line. It takes two arguments: path to fasta file and tne name of the new file. It returns new fasta file in the working directory.

 **Example:**

```
convert_multiline_fasta_to_oneline(path_to_fastа, name_of_ouput_file)
```

The `parse_blast_output` takes blast result file and creates new file. New file in the working directory contains the first names of genes in QUERY. It takes two arguments: path to fasta file and tne name of the new file. 

**Example:**

```
parse_blast_output(path_to_blast_result_file, name_of_ouput_file)
```