import sys


ALPHABET_FOR_DNA = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
ALPHABET_FOR_RNA = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
DICT_ALPHABET_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
DICT_ALPHABET_RNA = {'U': 'A', 'A': 'U', 'G': 'C', 'C': 'G',
                     'u': 'a', 'a': 'u', 'g': 'c', 'c': 'g'}
DICT_RNA_TO_PROTEIN = {  # dictionary triplet - amino acid
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L',
    'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I',
    'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S',
    'AGC': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y',
    'UAA': 'stop', 'UAG': 'stop', 'UGA': 'stop', 'CAU': 'H',
    'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'UGG': 'W',
    'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C', 'CGU': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'AUG': 'M'
    }


def transcribe(seqs: list[str]) -> list[str]:
    """
    Transcribe DNA-sequences into RNA-sequences

    :param seqs: list of arguments (str) DNA-sequences
    :return: list of arguments (str) DNA-sequences transcribed into RNA-sequences
    """
    output = []
    for seq in seqs:
        if not is_dna(seq):
            sys.exit(f'Sequence {seq} is not DNA')
        output.append(seq.replace('T', 'U').replace('t', 'u'))
    return output


def reverse(seqs: list[str]) -> list[str]:
    """
    Reverse DNA or RNA sequences

    :param seqs: list of arguments (str) DNA or RNA sequences
    :return: list of arguments (str) DNA or RNA sequences
    """
    output = []
    for seq in seqs:
        output.append(seq[::-1])
    return output


def is_dna(seq: str) -> bool:
    """
    Check if sequence is DNA

    :param seq: list of arguments (str) DNA or RNA sequences
    :return: - bool: True if the sequence is DNA sequence, False if it is RNA
    """
    unique_chars = set(seq)
    if unique_chars <= ALPHABET_FOR_DNA:
        return True
    else:
        return False


def complement(seqs: list[str]) -> list[str]:
    """
    Complement DNA or RNA sequences to RNA or DNA sequence

    :param seqs: list of arguments (str) DNA or RNA sequences
    :return: list of arguments (str) DNA or RNA sequences
    """
    output = []
    for seq in seqs:
        output_seq = []
        if is_dna(seq):  # if DNA
            for nucleotide in seq:
                output_seq.append(DICT_ALPHABET_DNA[nucleotide])
        else:  # if RNA
            for nucleotide in seq:
                output_seq.append(DICT_ALPHABET_RNA[nucleotide])
        output.append(''.join(output_seq))
    return output


def reverse_complement(seqs: list[str]) -> list[str]:
    """
    Reverse amd complement DNA or RNA sequences to RNA or DNA sequence

    :param seqs: list of arguments (str) DNA or RNA sequences
    :return: list of arguments (str) DNA or RNA sequences
    """
    output = (reverse(complement(seqs)))
    return output


def translation_rna_protein(seqs: list[str]) -> list[str]:
    """
    Translate RNA sequences protein

    :param seqs: list of arguments (str) RNA sequences
    :return: list of arguments (str) protein sequences
    """
    output = []
    for seq in seqs:
        output_process = []
        if len(seq) % 3 != 0 or ('T' in seq) or ('t' in seq):
            sys.exit(f'It is not possible to convert a given nucleic '
                     f'acid sequence {seq} into a protein.')
        else:
            triplets = [seq.upper()[i:i + 3] for i in range(0, len(seq), 3)]
            for triplet in triplets:
                output_process.append(DICT_RNA_TO_PROTEIN[triplet])
            output.append(''.join(output_process))
    return output


def translation_dna_protein(seqs: list[str]) -> list[str]:
    """
    Translate DNA sequences protein

    :param seqs: list of arguments (str) DNA sequences
    :return: list of arguments (str) protein sequences
    """
    output = []
    for seq in seqs:
        if len(seq) % 3 != 0 or 'U' in seq or 'u' in seq:
            output.append(f'To translate DNA into protein, please enter the '
                          f'DNA sequence, the sequence {seq} is not suitable.')
        else:
            output = translation_rna_protein(transcribe(seqs))
    return output


def is_nucleic_acid(seq: str) -> bool:
    """
    Checking for unrecorded characters

    :param seq: DNA or RNA sequences
    """
    unique_chars = set(seq)
    if unique_chars <= ALPHABET_FOR_DNA or unique_chars <= ALPHABET_FOR_RNA:
        return True
    else:
        sys.exit(f'I beg your pardon, the {seq} sequence is probably '
                 f'not nucleic acid \N{loudly crying face}')


def test_is_nucleic_acid(seqs: list[str]):
    """
    Checking nucleic acid if there are unrecorded characters

    :param seqs: list of arguments (str) DNA or RNA sequences
    """
    for seq in seqs:
        if is_nucleic_acid(seq):
            pass
