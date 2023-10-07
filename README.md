# Homework 5
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
  * [Usage example](#examples)

◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆

## A few words about the structure ✿ <a name="structure"></a> 
```
-/
 |- na_protein_fastq.py # (imports & 3 functions)
 |- README.md
 |- modules/
       |- dna_rna_tools.py
       |- protein_tools.py
       |- fastq_tool.py
```

There are 3 functions in the `na_protein_fastq.py` script: main function HW 3 (`run_dna_rna_tools.py`), main function HW 4 (`run_protein_tools.py`) and main function HW 5 (`run_fastq_tool.py`). The other functions are imported in the na_protein_fastq.py script, but are defined in the scripts in the `modules` folder.

### ᶘಠᴥಠᶅ What is in na_protein_fastq.py? ᶘಠᴥಠᶅ  <a name="what_is"></a> 
* The function `run_dna_rna_tools` takes as input an arbitrary number of arguments with DNA or RNA sequences (`str`) and the name of the procedure to be executed (always the last argument, `str`). After that the command performs the specified action on all submitted sequences using `modules/dna_rna_tools.py`. If one sequence is submitted, a string with the result is returned. If several sequences are passed - a list of strings is returned.
* The function `protein_tool.py` takes as input the name of the procedure and the sequence of amino acids, or two sequences, in the case of some procedures. This chair is designed for processing protein sequences using `modules/protein_tools.py`. The function returns a string that reflects the results of working with the sequence.
* Utility for working with fastq-sequences `run_fastq_tool.py` takes 4 arguments as input: `seqs`, `gc_bounds`, `length_bounds`, `quality_threshold`.
  1. `seqs` is a dictionary consisting of fastq sequences. The structure is as follows. Key - a string, the name of the sequence. Value - a tuple of two strings: sequence and quality. This is essentially the contents of the fastq file, but you and I haven't gone through reading files yet. So we'll use the python dictionary for now. Then it will be enough to add reading files and writing them to a dictionary of this kind, so that everything works from start to finish. The example_data.py script has an example for you to debug.
  2. `gc_bounds` - GC composition interval (in percent) for filtering (default is (0, 100), i.e. all sequences are saved). If a single number is passed in the argument, it is considered to be the upper bound. Examples: gc_bounds = (20, 80) - save only reids with GC composition from 20 to 80%, gc_bounds = 44.4 - save reids with GC composition less than 44.4%.
  3. `length_bounds` - length interval for filtering, everything is similar to gc_bounds, but by default it is equal to (0, 2**32).
  4. `quality_threshold` - threshold value of average quality of reid for filtering, default value is 0 (phred33 scale). Sequences with average quality of all nucleotides below the threshold are discarded.


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
**3. fastq-seqbences selection**
``` python
seqs = {
    '@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804:1:SRR292678:1:1101:47176:47176': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB')
}
print(run_fastq_tool(seqs, gc_bounds=100, length_bounds=(0, 2 ** 32), quality_threshold=35))
```

◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
