import sys
import modules.dna_rna_tools as dna_rna_tools
import modules.protein_tools as protein_tools
import modules.fastq_tool as fastq_tool


def run_dna_rna_tools(*function_input: str) -> list[str] or str:
    """
    Main function that is used to get protein sequence(s) and command. It performs a given action

    :param function_input: arbitrary number of arguments with DNA or RNA sequences (str)
    and the name of the procedure to be executed (last argument, str)
    :return: If one sequence is submitted, a string with the result is returned.
    If several sequences are passed - a list of strings is returned.
    """
    *seqs, function = function_input

    dna_rna_tools.test_is_nucleic_acid(seqs)  # Nucleic acid test

    dict_of_functions = {'transcribe': dna_rna_tools.transcribe,
                         'reverse': dna_rna_tools.reverse,
                         'reverse_complement': dna_rna_tools.reverse_complement,
                         'complement': dna_rna_tools.complement,
                         'translation_rna_protein': dna_rna_tools.translation_rna_protein,
                         'translation_dna_protein': dna_rna_tools.translation_dna_protein,
                         }
    answer = dict_of_functions[function](seqs)

    if len(answer) == 1:
        return answer[0]
    else:
        return answer


def run_protein_tools(*args: str) -> str:
    """
    Main function that is used to get protein sequence(s) and command. It performs a given action

    Arguments:
    - args (str): amino acid sequence(s) and command.
    The input must use the single letter amino acid code
    The last element of the string must be the command

    Returns:
    - result (str): the result of a given sequence processing
    """
    *sequences, action = args
    sequences = [seq.upper() for seq in sequences]
    commands = {'calculate_amino_acid_percentages': protein_tools.calculate_amino_acid_percentages,
                'classify_amino_acid': protein_tools.classify_amino_acid,
                'find_amino_acid_indices': protein_tools.find_amino_acid_indices,
                'counting_point_mutations': protein_tools.counting_point_mutations,
                'counting_molecular_weight': protein_tools.counting_molecular_weight,
                'get_occurrences': protein_tools.get_occurrences,
                'count_variant_rna': protein_tools.count_variant_rna,
                'determine_total_protein_charge': protein_tools.determine_total_protein_charge,
                'calculate_pI': protein_tools.calculate_pi}
    command = commands[action]
    for seq in sequences:
        if not protein_tools.is_protein(seq):
            print('Sequence is not protein', file=sys.stderr)
            sys.exit(1)
    return str(command(sequences[0], sequences[1])) if len(sequences) == 2 else str(command(sequences[0]))


def run_fastq_tool(seqs: dict, gc_bounds: tuple or int = (0, 100), length_bounds: tuple or int = (0, 2 ** 32),
                   quality_threshold: int = 0) -> dict:
    """
    Main function that is used to get protein sequence(s) and command. It performs a given action

    :param seqs: seqs is a dictionary consisting of fastq sequences
    :param gc_bounds: GC composition interval (in percent) for filtering (default is (0, 100)
    :param length_bounds: length interval for filtering (default is (0, 2 ** 32)
    :param quality_threshold: threshold value of average quality of reid for filtering, default value is 0
    :return: dictionary consisting of only those sequences that passed all conditions
    """
    out_dict = {}
    for seq in seqs:
        filter_output = (
                fastq_tool.count_gc(seqs[seq], gc_bounds) and
                fastq_tool.count_length(seqs[seq], length_bounds) and
                fastq_tool.count_quality(seqs[seq], quality_threshold)
        )
        if filter_output:
            out_dict[seq] = seqs[seq]
    return out_dict
