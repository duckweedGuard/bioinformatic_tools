def count_gc(seq: tuple, gc_bounds: tuple or int = (0, 100)) -> bool:
    """
    Function calculates GC composition and if value is in GC composition interval for filtering, it returns True,
    else False.

    :param seq: tuple of two strings: sequence and quality
    :param gc_bounds: GC composition interval (in percent) for filtering (default is (0, 100)
    :return: True if the obtained value is included in gc_bounds, else False
    """
    if type(gc_bounds) != tuple:
        gc_bounds = tuple([0, gc_bounds])
    gc_result = (seq[0].count('C') + seq[0].count('G')) / len(seq[0]) * 100
    if gc_result < gc_bounds[0]:
        return False
    elif gc_result > gc_bounds[1]:
        return False
    else:
        return True


def count_length(seq: tuple, length_bounds: tuple or int = (0, 2 ** 32)) -> bool:
    """
    Finds the length and if the obtained value is not included in the length interval for filtering, function
    returns False, else it returns True.

    :param seq: tuple of two strings: sequence and quality
    :param length_bounds: length interval for filtering, default (0, 2**32)
    :return: True if the obtained value is included in count_length, else False
    """
    if type(length_bounds) != tuple:
        length_bounds = tuple([0, length_bounds])
    length_result = len(seq[0])
    if length_result < length_bounds[0]:
        return False
    elif length_result > length_bounds[1]:
        return False
    else:
        return True


def count_quality(seq: tuple, quality_threshold: int = 0) -> bool:
    """
    Finds the average quality of the reed

    :param seq: tuple of two strings: sequence and quality
    :param quality_threshold: threshold value of average quality of reid for filtering, default value is 0
    :return: True if average quality less than quality_threshold, else False
    """
    quality_result = 0
    for sign in seq[1]:
        quality_result += (ord(sign) - 33)
    quality_out = quality_result / len(seq[1])
    if quality_out < quality_threshold:
        return False
    else:
        return True


def read_fastq(input_path: str, input_filename: str) -> dict:
    """
    Function to read a file by a given path and file name. Then the data is translated into a dictionary

    :param input_path: given path where input file is
    :param input_filename: input file name
    :return: dictionary consisting of fastq sequences
    """
    if not os.path.isfile(os.path.join(input_path, input_filename)):
        raise ValueError('No such fastq file in directory')
    with open(os.path.join(input_path, 'fastq.fastq'), mode='r') as file:
        fastq_dict = {}
        for line in file:
            if ' ' in line and '+' not in line:
                seq_line = file.readline().strip()
                file.readline()
                quality_line = file.readline().strip()
                fastq_dict[line] = tuple([seq_line, quality_line])
        return fastq_dict


def write_output_fastq(out_dict: dict, input_filename: str,
                       output_data_dir: str = 'fastq_filtrator_resuls',
                       output_filename: str = None):
    """
    Function for writing filtered fastq sequences

    :param out_dict: post-work filtered dictionary
    :param input_filename: input file name
    :param output_data_dir: path where the file with filtered sequences will be saved
    :param output_filename: name file with filtered sequences will be saved
    """
    if not os.path.isdir(output_data_dir):
        os.mkdir(output_data_dir)
    if output_filename is None:
        output_filename = 'filtered_' + input_filename
    with open(os.path.join(output_data_dir, output_filename), mode='w') as file:
        for key in out_dict:
            file.write(key)
            file.write(out_dict[key][0] + "\n")
            file.write("+" + key)
            file.write(out_dict[key][1] + "\n")

