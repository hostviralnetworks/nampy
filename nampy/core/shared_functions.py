def test_kwarg(the_arg_string, the_kwarg_dict, the_allowed_values):
    """ A simple function to help set default values.

    Arguments:
     the_test_string: a desired key to check for in
                      the_kwarg_dict
     the_kwarg_dict: the dict that have a preferred value
     the_allowed_values: a list of allowed values.
                         the first entry is the default.
     

    """
    if the_arg_string in the_kwarg_dict.keys():
        if the_kwarg_dict[the_arg_string] in the_allowed_values:
            the_set_value = the_kwarg_dict[the_arg_string]
        else:
            the_set_value = the_allowed_values[0]
    else:
        the_set_value = the_allowed_values[0]

    return the_set_value
