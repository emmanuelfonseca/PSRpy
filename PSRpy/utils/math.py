def represents_an_int(input_string):
    """
    Checks if string can be represented as a Python integer and returns a boolean.
    This function is used to evaluate if 

    taken from: 
        https://stackoverflow.com/questions/1265665/
        how-can-i-check-if-a-string-represents-an-int-without-using-try-except
    """

    is_an_integer = False

    try:
        value = int(input_string)
        is_an_interger = True

    except ValueError:
        pass

    return is_an_integer
