# Homework 6
## bioinformatic_tools ･ﾟ✧(=✪ ᆺ ✪=)･ﾟ✧

**ᶘ⊙ᴥ⊙ᶅ** Author: Orlova Victoria.

<img src="https://www.meme-arsenal.com/memes/6e7a90e11e31bbe40c15cdff7e442c92.jpg" width="200" height="200">

◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
## Table of contents

  * [A few words about the structure](#structure)
  * [What is in na_protein_fastq.py script?](#what_is)
  * [Module dna_rna_tools.py functional content](#module_dna_rna) 
  * [Module protein_tools.py functional content](#module_protein)
  * [Module fastq_tool.py functional content](#module_fastq)
  * [What is in bio_files_processor.py?](#what_is_bio)
  * [Usage example](#examples)

◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆

## A few words about the structure ✿ <a name="structure"></a> 
```
-/
 |- na_protein_fastq.py # (imports & 3 functions)
 |- bio_files_processor.py # (4 functions)
 |- README.md
 |- fasta.fasta
 |- modules/
       |- dna_rna_tools.py
       |- protein_tools.py
       |- fastq_tool.py
 |- data/
       |- short_example_gbk.gbk
       |- fastq.fastq
       |- start_pos_fasta.fasta
       |- example_blast_results.txt

```

There are 3 functions in the `na_protein_fastq.py` script: main function HW 3 (`run_dna_rna_tools.py`), main function HW 4 (`run_protein_tools.py`) and :sparkles:updated:sparkles: in HW6 main function HW5 (`run_fastq_tool.py`). The other functions are imported in the na_protein_fastq.py script, but are defined in the scripts in the `modules` folder.

:sparkles: New :sparkles: There are 4 functions in the `bio_files_processor.py` script: function `convert_multiline_fasta_to_oneline`, function `select_genes_from_gbk_to_fasta`, function `change_fasta_start_pos` and function `parse_blast_output`. Functions use file-based input and output. All input files lie in the `date` folder, output files are created in the root of the repository. Read the descriptions of functions for more details.


### ᶘಠᴥಠᶅ What is in na_protein_fastq.py? ᶘಠᴥಠᶅ  <a name="what_is"></a> 
* The function `run_dna_rna_tools` takes as input an arbitrary number of arguments with DNA or RNA sequences (`str`) and the name of the procedure to be executed (always the last argument, `str`). After that the command performs the specified action on all submitted sequences using `modules/dna_rna_tools.py`. If one sequence is submitted, a string with the result is returned. If several sequences are passed - a list of strings is returned.
* The function `protein_tool.py` takes as input the name of the procedure and the sequence of amino acids, or two sequences, in the case of some procedures. This chair is designed for processing protein sequences using `modules/protein_tools.py`. The function returns a string that reflects the results of working with the sequence.
* :sparkles: New :sparkles: Utility for working with fastq-sequences `run_fastq_tool.py` takes 7 arguments as input: `input_path`, `input_filename`, `gc_bounds`, `length_bounds`, `quality_threshold`,  `output_filename`, `output_data_dir`.
  1. `input_path` - given path where input file is (in our example `input_path = 'data'`) 
  2. `input_filename` - given input file name
  3. `gc_bounds` - GC composition interval (in percent) for filtering (default is (0, 100), i.e. all sequences are saved). If a single number is passed in the argument, it is considered to be the upper bound. Examples: gc_bounds = (20, 80) - save only reids with GC composition from 20 to 80%, gc_bounds = 44.4 - save reids with GC composition less than 44.4%.
  4. `length_bounds` - length interval for filtering, everything is similar to gc_bounds, but by default it is equal to (0, 2**32).
  5. `quality_threshold` - threshold value of average quality of reid for filtering, default value is 0 (phred33 scale). Sequences with average quality of all nucleotides below the threshold are discarded.
  6. `output_filename` - name file with filtered sequences will be saved
  7. `output_data_dir` - path where the file with filtered sequences will be saved, if `None`, will create `fastq_filtrator_resuls` folder


**(∿°○°)∿ .・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭**


## ʕ•ᴥ•ʔﾉ♡ Module dna_rna_tools.py functional content   ʕ•ᴥ•ʔﾉ♡  <a name="module_dna_rna"></a> 
### Functions and brief description of how they work: ✿

* The function `test_is_nucleic_acid` takes as input the list of arguments (str) received at the input of the function run_dna_rna_tools. The command then checks whether the input sequences are nucleic acid using the is_nucleic_acid function when checking the elements of the list (sequences). If the sequence is not DNA or RNA, the program is aborted with a comment returned.
* `is_nucleic_acid` is a function to check a sequence for unrecorded characters. It takes a string - sequence - as input. After that it is processed for cases when the sequence cannot exist because it contains T and U at the same time. Returns True if the sequence is a nucleic acid and a comment if the sequence is not a nucleic acid.
* The `transcribe` function takes as input the list of arguments (str) received at the input of the function run_dna_rna_tools. After operation, it returns DNA-sequences transcribed into RNA-sequences.
* The `reverse` function takes as input the list of arguments (str) received at the input of the function run_dna_rna_tools. After operation, it returns reversed DNA-sequences or RNA-sequences.
* The `complement` function takes as input the list of arguments (str) received at the input of the function run_dna_rna_tools. After running, it returns the complemented sequences. First, for a sequence in the list, the DNA or RNA membership is checked using the is_dna function. DNA sequence nucleotides are replaced with complementary ones using the DICT_ALPHABET_DNA dictionary for DNA nucleotides, and RNA sequence nucleotides are replaced with complementary ones using the DICT_ALPHABET_RNA dictionary for RNA nucleotides.
* `is_dna` is a function to determine if the input sequence is DNA. It takes a sequence (str) as input and returns True if the sequence is a DNA sequence or returns False if the sequence is an RNA sequence.
* The `reverse_complement` function takes as input the list of arguments (str) received at the input of the function run_dna_rna_tools. It returns the reverse complementary sequence. First the complementary one is created by running the transcribe function, then the sequence is reversed by the reverse function.
* The `translate_rna_protein` function takes as input the list of arguments (str) received at the input of the function run_dna_rna_tools. It translates the RNA sequence into a protein. The DICT_RNA_TO_PROTEIN dictionary is used, where amino acids (one-letter code) correspond to RNA codons. First, a check is made that the sequence is an RNA sequence and that it can be divided into codons. Next, the amino acid corresponding to the codon is found and added to the output_str string, after all triplets have been processed, the output_str is added to the output list.
* The `translation_dna_protein` function takes as input the list of arguments (str) received at the input from the `run_dna_rna_tools`. It performs translation of a DNA sequence into a protein. The function first checks whether the sequence is DNA and whether it can be divided into codons. Then the DNA sequence is transcribed into an RNA sequence by the function transcribe, after which the RNA sequence is translated into a protein using the function translate_rna_protein.


## ʕ•ᴥ•ʔﾉ♡ Module protein_tools.py functional content   ʕ•ᴥ•ʔﾉ♡ <a name="module_protein"></a> 
*Some rules:*
- The sequence should be composed of single-letter names of amino acids
- The sequence can be passed in uppercase or lowercase
- The sequence must contain proteinogenic amino acids, without their modifications

### Functions and brief description of how they work: ✿
* The `count_variant_rna` function accepts a protein sequence (str) as input. After that, the function counts the number of possible RNA variants that can be a matrix for the synthesis of a given amino acid sequence. The result is returned, the number of possible RNAs (int)
* The `determine_total_protein_charge` function accepts a protein sequence (str) as input. After that, the function determines whether a given amino acid sequence is positively negatively charged or not charged. The result is returned, the string `negative`, `positive`, `neutral`
* The `calculate_pi` function accepts a protein sequence (str) as input. After that, the function calculates the approximate value of the isoelectric point (pI) of a given amino acid sequence. The result is returned, the isoelectric point (float)
* The `counting_point_mutations` function takes two protein sequences (str) as input. Input sequences must have the same length. After that, the function counts the number of mutations - amino acid substitutions, the result returned is the number of mutations (int).
* The `counting_molecular_weight` function takes a protein sequence (str) as input. The function then counts the molecular weight of the input protein sequence, the result returned is the number of molecular weight (int).
* The function `get_occurrences` takes as input two protein sequences (str) - seq1 and seq2. After that, the command counts the number of occurrences without intersection of the second protein sequence into the first one. The result is returned as a string (str), where across the space: the number of occurrences of the second string into the first string, followed by the occurrence indices (1 and beyond) across the space.
* The `calculate_amino_acid_percentages` function takes a protein sequence (str) as input. The command then calculates the percentage of amino acids in the protein, and finally returns the result as a string (str), where the amino acid and its percentage in the protein are written with a colon. The function rounds the result to two decimal places and sorts it in descending order.
* The function `classify_amino_acid` takes a protein sequence (str) as input. The command returns the calculated percentage of non_charged, acidic and basic amino acids in the protein as a string (str). The function rounds the result to two decimal places.
* The `find_amino_acid_indices` function takes as input a protein sequence and an amino acid (str) - seq and amino_acid. After that the command searches for indices of amino acid occurrence in protein and ultimately returns the result as a string (str) with all the indices.
Using this function you may encounter an error **ValueError: Amino acid not found**. This error occurs if the entered amino acid is not contained in the protein sequence or it is written in lower case.


## ʕ•ᴥ•ʔﾉ♡ Module fastq_tool.py functional content   ʕ•ᴥ•ʔﾉ♡ <a name="module_fastq"></a> 
### Functions and brief description of how they work: ✿
* The `count_gc` function takes a tuple of two strings: sequence and quality and gc_bounds as input. Function calculates GC composition and if the obtained value is not within the GC composition interval for filtering, it returns False, else function returns True.
* The `filter_length` function takes a tuple of two strings: sequence and quality and length_bounds as input. Finds the length and if the obtained value is not included in the length interval for filtering, function returns False, else it returns True.
* The `filter_quality` function takes a tuple of two strings: sequence and quality and quality_threshold as input. Finds the average quality of the reed and if the obtained value is greater than the threshold value of the average quality of the reed for filtering, function returns False, else it returns True.
* :sparkles: New :sparkles: The `read_fastq` function takes a input_path (str) and input_filename (str) as input. Read file and translate data into a dictionary consisting of fastq sequences. The structure is as follows. Key - a string, the name of the sequence. Value - a tuple of two strings: sequence and quality.
* :sparkles: New :sparkles: The `write_output_fastq` function takes a out_dict (post-work filtered dictionary), input_filename (str), output_data_dir (str), output_filename (str). Function writes filtered fastq seqbenses into a fastq-file.


### ᶘಠᴥಠᶅ What is in bio_files_processor.py? ᶘಠᴥಠᶅ  <a name="what_is_bio"></a> :sparkles: New :sparkles:
There are 4 functions in the `bio_files_processor.py` script: function `convert_multiline_fasta_to_oneline`, function `select_genes_from_gbk_to_fasta`, function `change_fasta_start_pos` and function `parse_blast_output`. 
* The `convert_multiline_fasta_to_oneline` takes an input_fasta (str) as imput file name. Reads an input fasta file in which the sequence (DNA/RNA/protein/ ... ) can be split into several lines, then saves it into a new fasta file in which each sequence fits into one line. The `.fasta` extension is appended to the `output_fasta` name, if `output_fasta = None`, then 'converted_' is appended to the input_fasta file name.
* The 'select_genes_from_gbk_to_fasta' takes an input_gbk (str), genes (list of str), n_before (int, default = 1), n_after (int, default = 1), output_fasta (default = 'output_neighboor_genes_fasta'). Function select n_before number of genes before and n_after number of genes after each gene of interest from .gbk file and save their protein sequence (translation) to a fasta file. Use `bring_indexes` function below
* `bring_indexes` function receives as input a list of all genes found in the file (all_genes), a list of target genes - genes of interest (matching_genes), list for answer (out_neighboor_genes), number of genes before (>0, default value is 1), number of genes after (>0, default value is 1). Then selects genes neighboring the genes of interest, located at a given distance before or after the gene of interest. In case the gene of interest is closer to the beginning/end than the number of neighbors we are looking for, it selects only the possible ones. If two genes of interest are next/ close to each other, they are not included in the final out_neighboor_genes list.
* The `change_fasta_start_pos` function takes an input_fasta (str), shift (int), output_fasta (str) as input. Function shift the initial position in the one-line fasta with one record. If `abs(shift)` is bigger than fasta-record sequense, then there is an correction by division.
* The `parse_blast_output` function takes an input_file (str), output_file (str) (if `output_file = None`, then `output_file = 'parse_blast_blast_results.txt'`) as input. Function reads the given txt file, for each QUERY select the first line from the Description column. Save the list of obtained proteins to a new file in a single column sorted alphabetically.



◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
### Usage example: <a name="examples"></a> 


**☆ ☆ ☆ :. o(≧▽≦)o .: ☆ ☆ ☆**

**1. Operations with nucleic acids sequences** 
``` python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('ATG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
run_dna_rna_tools('Augccg', 'AauAuGCCaCaCaC', 'translation_rna_protein') # ['MP', 'NMPH']
run_dna_rna_tools('ATg', 'AaTttGCCaCaC', 'translation_dna_protein') # ['M', 'NLPH']
```

**☆ ☆ ☆ :. o(≧▽≦)o .: ☆ ☆ ☆**

**2. Operations with protein sequences**
``` python
protein_tool('TATAQQQWRVVTDDDA', 'count_variant_rna') # '25165824'
protein_tool('TDDDTEQQWRVVTDDDA', 'determine_total_protein_charge') # 'negative'
protein_tool('TKKKKTDDDA', 'calculate_pi') # '7.225555555555555'
run_protein_tools('ASQG', 'AMQR', 'counting_point_mutations') # 2
run_protein_tools('ASQGAMQR', 'counting_molecular_weight') # 847
run_protein_tools('ASQRGARWQRMQR', 'QR', 'get_occurrences') # 'Number of occurrences: 3; indexes: 3, 9, 12'
run_protein_tools('ADNNDQD', 'calculate_amino_acid_percentages') # 'D: 42.86, N: 28.57, A: 14.29, Q: 14.29'
run_protein_tools('ARNDCQ', 'classify_amino_acid') # 'non_charged: 66.67, acidic: 16.67, basic: 16.67'
run_protein_tools('ARNDCQA', 'A', 'find_amino_acid_indices') # '1, 7'
```

**☆ ☆ ☆ :. o(≧▽≦)o .: ☆ ☆ ☆**

 **3. fastq-seqbences selection** :sparkles: New :sparkles:
``` python
print(run_fastq_tool('data', 'fastq.fastq', 40, 60, 20))
```
After the function is run, the `fastq.fastq` file from the folder`data` will be read, then filtering will be done and the answer will be written to the `fastq_filtrator_resuls` folder (if it does not exist, it will be created) under the entered  `output_filename` (if the `output_filename` was not entered, it will be named as follows: `filtered_` + input_filename, in example `filtered_fastq.fastq`)


**☆ ☆ ☆ :. o(≧▽≦)o .: ☆ ☆ ☆**

 **4. Operations with bio_files_processor.py** :sparkles: New :sparkles:
 
**4.1. convert_multiline_fasta_to_oneline**

``` python
convert_multiline_fasta_to_oneline('fasta.fasta')
```


**4.2. select_genes_from_gbk_to_fasta**

``` python
nput_gbk = 'data/short_example_gbk.gbk'
genes = ['mngR', 'pxpB']
n_before = 2
n_after = 2
select_genes_from_gbk_to_fasta(input_gbk, genes, n_before, n_after)
```

**4.3. change_fasta_start_pos**

``` python
input_fasta = 'data/start_pos_fasta.fasta'
shift = -1
output_fasta = 'changed_start_pos_fasta.fasta'
change_fasta_start_pos(input_fasta, shift, output_fasta)
```


**4.4. parse_blast_output**
``` python
input_file = 'data/example_blast_result.txt'
parse_blast_output(input_file)
```



◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
