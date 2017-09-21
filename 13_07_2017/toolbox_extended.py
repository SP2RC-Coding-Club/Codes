import numpy as np

# Check if a number is within some interval of any other value in a list.
# If it is, then it gets discarded. If it's not, it gets appended to the list.
# Example:
# If you want to check if 3 is within 2 of any element in the list [5,4,10],
# you would input check(3, [5,4,10], 2) and obtain the list [5,4,10] as output
# since 3 is within +/-2 of the element 4.
# If you want to check if 3 is within 0.1 of any element in the list [5,4,10],
# you input check(3, [5,4,10], 0.1) and obtain the list [5,4,10,3] as output
# since 3 is not within +/-0.1 of any element in the list [5,4,10].
def check(value, value_list, difference):
    n = True
    if len(value_list) == 0:
        value_list.append(value)
    else:
        for x in value_list:
            if np.abs(np.real(x) - np.real(value)) < difference and \
               np.abs(np.imag(x) - np.imag(value)) < difference:
                   n = False
            else:
                pass
        if n == True:
            value_list.append(value)
    return value_list

# This function converts a list of lists into a numpy array. It only takes the
# list of lists as input, and returns the array as output. If the lists inside 
# the list are of unequal lengths, it fills up the lines with None so that all 
# lines in the output array are of equal length.
# Example input:
# a = [[1,3,4], [2,1], [2,3,4,7]]
# Output:
# array([[1, 3, 4, None],
#        [2, 1, None, None],
#        [2, 3, 4, 7]], dtype=object)
def list_to_array(the_list):
    length = len(sorted(the_list, key=len, reverse=True)[0])
    array = np.array([xi+[None]*(length-len(xi)) for xi in the_list])
    array = array.astype(np.complex)
    return array