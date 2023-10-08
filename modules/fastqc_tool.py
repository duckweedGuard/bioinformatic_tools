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
    if quality_out > quality_threshold:
        return False
    else:
        return True
