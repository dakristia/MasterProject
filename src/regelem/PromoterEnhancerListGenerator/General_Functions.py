def check_in_range(value_to_check1, value_to_check2, start_range1, end_range1, start_range2, end_range2):
    end_range1 += 1; end_range2 += 1
    return value_to_check1 in range(start_range1,end_range1) or value_to_check1 in range(start_range2,end_range2) or value_to_check2 in range(start_range1,end_range1) or value_to_check2 in range(start_range2,end_range2)

def nums_to_str_range(rangeFrom,rangeTo):
    return str(rangeFrom) + "-" + str(rangeTo)

def get_digit(number, n):
    return number // 10**n % 10

def round_up_and_down(num : int, interval : int):
    """_summary_

    Args:
        num (_int_): _description_
        interval (_int_): _description_

    Returns:
        tuple: tuple with number rounded down and up (in that order)
    """
    #Left part of tuple is rounded down, right part is rounded up
    rest = num % interval
    minimum = num - rest

    maximum = minimum + interval
    return (minimum,maximum)