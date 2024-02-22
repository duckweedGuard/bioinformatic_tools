# Homework 14
## bioinformatic_tools ･ﾟ✧(=✪ ᆺ ✪=)･ﾟ✧

**ᶘ⊙ᴥ⊙ᶅ** Author: Orlova Victoria.

<img src="https://www.meme-arsenal.com/memes/6e7a90e11e31bbe40c15cdff7e442c92.jpg" width="200" height="200">

◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
## Table of contents

  * [A few words about the structure](#structure)
  * [filter_fastq.py](#what_is)
  * [BiologicalSequencesWorld.py](#what_is_bio)
  * [Usage example](#examples)

◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆

## A few words about the structure ✿ <a name="structure"></a> 
```
-/
 |- filter_fastq.py
 |- README.md
 |- BiologicalSequencesWorld.py
 |- bio_files_processor.py # just ignore it
 |- data/
       |- fastq.fastq

```


### ᶘಠᴥಠᶅ filter_fastq.py? ᶘಠᴥಠᶅ  <a name="what_is"></a> 
* Utility for working with fastq-sequences `filter_fastq.py` takes 7 arguments as input: `input_path`, `input_filename`, `gc_bounds`, `length_bounds`, `quality_threshold`,  `output_filename`, `output_data_dir`.
  1. `input_path` - given path where input file is. 
  2. `input_filename` - given input file name
  3. `gc_bounds` - GC composition interval (in percent) for filtering (default is (0, 100), i.e. all sequences are saved). If a single number is passed in the argument, it is considered to be the upper bound. Examples: gc_bounds = (20, 80) - save only reids with GC composition from 20 to 80%, gc_bounds = 44.4 - save reids with GC composition less than 44.4%.
  4. `length_bounds` - length interval for filtering, everything is similar to gc_bounds, but by default it is equal to (0, 2**32).
  5. `quality_threshold` - threshold value of average quality of reid for filtering, default value is 0 (phred33 scale). Sequences with average quality of all nucleotides below the threshold are discarded.
  6. `output_filename` - name file with filtered sequences will be saved
  7. `output_data_dir` - path where the file with filtered sequences will be saved, if `None`, will create `fastq_filtrator_resuls` folder


**(∿°○°)∿ .・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭・.・。.・゜✭**


### ᶘಠᴥಠᶅ BiologicalSequencesWorld? ᶘಠᴥಠᶅ  <a name="what_is_bio"></a> 
An abstract BiologicalSequence class that defines the following interface:

*    Working with the len function
*    Ability to retrieve elements by index and make sequence slices (similar to strings)
*    Print output in a convenient form and the ability to convert to string
*    Ability to check the sequence alphabet for correctness

The NucleicAcidSequence class:

*    This class implements the BiologicalSequence interface
*    This class has a new complement method that returns the complementary sequence
*    This class has a new method gc_content that returns the GC composition (no difference whether it is a percentage or a fraction)

The successor classes to NucleicAcidSequence: DNASequence and RNASequence

*    DNASequence must have a transcribe method that returns the transcribed RNA sequence
*    These classes should not have public complement methods and a method to check the alphabet, as these should already be implemented in NucleicAcidSequence.

An AminoAcidSequence class:

*    This class implements the BiologicalSequence interface
*    Method for calculating the molecular weight of a protein.



◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
### Usage example: <a name="examples"></a> 


**☆ ☆ ☆ :. o(≧▽≦)o .: ☆ ☆ ☆**

**1. FastQ-filter using Biopython** 
``` python
input_path = "data/fastq.fastq"
output_filename = 'filtered_fastq.fastq'
filter_fastq(input_path, gc_bounds = (45, 100), length_bounds = (0, 2 ** 32), quality_threshold = 0)
```

**☆ ☆ ☆ :. o(≧▽≦)o .: ☆ ☆ ☆**

**2. Operations with protein sequences**

``` python
from BiologicalSequencesWorld import *
my_dna = BiologicalSequence('AGTCC')
my_dna.type_definition()
my_dna
# The sequence is: AGTCC, type is DNA
```

``` python
my_dna = NucleicAcidSequence('ATTGCAT')
my_dna.complement()
# The sequence is: TAACGTA, type is DNA
```

``` python
my_dna.gc_content()
# 28.571
```

``` python
my_dna = DNASequence('AGTCC')
my_dna.transcribe()
# The sequence is: AGUCC, type is RNA
```

``` python
my_rna = RNASequence('AUUGCAU')
my_rna.transcribe_reverse()
The sequence is: TACGTTA, type is DNA
```

``` python
my_protein = AminoAcidSequence('FASIT')
my_protein.counting_molecular_weight()
# 537
```


◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆◇◆
