from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os


def filter_fastq(input_path: str, gc_bounds: tuple or int = (0, 100),
                 length_bounds: tuple or int = (0, 2 ** 32),
                 quality_threshold: int = 0,
                 output_filename=None,
                 output_data_dir: str = 'filter_fastq_results'):
    """
    Main function that is used to get fastq sequence(s) and command. It performs a given action
    :param input_path: given path where input file is
    :param gc_bounds: GC composition interval (in percent) for filtering (default is (0, 100)
    :param length_bounds: length interval for filtering (default is (0, 2 ** 32)
    :param quality_threshold: threshold value of average quality of reid for filtering, default value is 0
    :param output_filename: name file with filtered sequences will be saved
    :param output_data_dir: path where the file with filtered sequences will be saved
    """
    if not os.path.isdir(output_data_dir):
        os.mkdir(output_data_dir)
    if output_filename is None:
        output_filename = 'filtered_fastq'

    with open(input_path) as handle, open(os.path.join(output_data_dir, output_filename), mode='w') as file:
        record = SeqIO.parse(handle, "fastq")
        for lin in record:
            out = 0
            seq = lin.seq

            if type(length_bounds) != tuple:
                length_bounds = tuple([0, length_bounds])
            length_result = len(lin.letter_annotations['phred_quality'])
            if not length_bounds[0] <= length_result <= length_bounds[1]:
                out = 1

            if not isinstance(gc_bounds, tuple):
                gc_bounds = (0, gc_bounds)
            gc_result = gc_fraction(seq) * 100
            if not gc_bounds[0] <= gc_result <= gc_bounds[1]:
                out = 1

            score = lin.letter_annotations['phred_quality']
            score = sum(score) / len(seq)
            if not score >= quality_threshold:
                out = 1

            if out == 0:
                file.write(lin.format("fastq"))
