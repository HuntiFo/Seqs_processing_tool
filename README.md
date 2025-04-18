# Seqs_processing_tool

This is a bioinformatic package for processing nucleotide sequences, which consists three classes to work with biological sequences and function `filter_fastq`.

## Installation

To install Seqs_processing_tool you need to download the repository from GitHub

```
git clone https://github.com/HuntiFo/Seqs_processing_tool.git
cd Seqs_processing_tool
```

New directory Seqs_processing_tool contains:
1) the main script `main.py` (three classes to work with biological sequences and fucntion `filter_fastq`).
2) bio_files_processor.py 
and README.md.

## Running instructions

### Main function

To work with biological sequences create an object belonging to one of three classes: DNASequence, RNASequence, AminoAcidSequence with an input argument sequence as a string. Each class contains function `check_alphabet` to define input sequence.

1) Classes DNASequence and RNASequence both contain methods: `reverse`(return reverse sequence), `complement`(return complement sequence) and `reverse_complement`(return reverse-complement sequence). Class DNASequence has unique function `transcribe` to convert DNA into RNA.

2) Class AminoAcidSequence has method `calculate_weight` to calculate in weight of protein in g/mol.

**Example:**
```
dna = DNASequence('ACTGC')
dna.complement()
dna.transcribe()
```

B) `filter_fastq`
To run `filter_fastq` open the main.py, write the path to the fastq file, the name to the new file and next additional options. Call the function. 
The function have five arguments: input_fastq,output_fastq, gc_bounds, length_bounds, quality_threshold.
1) `input_fastq` - str, the path to fastq reads.
2) `output_fastq` - str, the name of a new file with selected fastq reads.
3) `gc_bounds` - str or tuple, the GC interval of the composition (in percent) for filtering (by default (0, 100)). If you pass a single number to the argument, it is assumed that this is the upper bound. Reads with GC composition higher than upper bound and lower than lower bound are discarded.
3) `length_bounds` - the length interval for filtering (by default it is (0, 2**32)). Reads with length longer than upper bound and shorter than lower bound are discarded.
4) `quality_threshold` - the threshold value of the average read quality for filtering (0 by default) (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded.
All the described intervals include both upper and lower bounds.

 **Example:**
```
filter_fastq('./dir1/read.fastqc', './dir1/read_filtrated.fastqc', gc_bounds = (10,90), length_bounds = (30,90), quality_threshold = 5)
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