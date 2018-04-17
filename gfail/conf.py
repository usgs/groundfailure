#!/usr/bin/env python

# stdlib imports
import os.path

# third party libraries
from configobj import ConfigObj
from validate import Validator, VdtTypeError


def __getCustomValidator():
    '''
    Return a Validator object with the custom types we have defined here.

    Returns:
        Validator object with custom types embedded.
    '''
    fdict = {
        'file_type': __file_type,
        'path_type': __path_type,
    }

    validator = Validator(fdict)
    return validator


def __file_type(value):
    '''
    Describes a file_type from the ShakeMap config spec.
    A file_type object is simply a string that must be a valid file on the
    system.

    Args:
        value (str): Path to a file on the local system.

    Returns:
        str: Input string, if a valid file name.
    '''
    if not os.path.isfile(value):
        raise VdtTypeError(value)
    return value


def __path_type(value):
    '''
    Describes a path_type from the groundfailure config spec.
    A path_type object is simply a string that must be a valid file OR
    directory on the system.

    Args:
        value (str): Path to a file or directory on the local system.

    Returns:
        str: Input string, if a valid file/directory name.
    '''
    if not os.path.isfile(value) and not os.path.isdir(value):
        raise VdtTypeError(value)
    return value


def __filterResults(result):
    # TODO: this function has a problem where some error messages are
    # duplicated...?
    errormsg = ''
    for key, value in result.items():
        if isinstance(value, dict):
            tmpmsg = __filterResults(value)
            errormsg += tmpmsg
        else:
            if not isinstance(value, bool):
                errormsg += ("Parameter %s failed with error '%s'\n"
                             % (key, value.args))
            else:
                if not value:
                    errormsg += ("Parameter %s was not specified correctly.\n"
                                 % (key))

    return errormsg


def correct_config_filepaths(input_path, config):
    """
    Takes an input filepath name and pre-pends it to all file locations within
    the config file. Individual locations are put into the config.  Don't have
    to put entire filepath location for each layer. Works by looping over
    config dictionary and subdictionary to fine locations named 'file'.

    Args:
        input_path (str): Path that needs to be appended to the front of all
            the file names/paths in config.
        config (ConfigObj): Object defining the model and its inputs.

    Returns:
        config dictionary with complete file paths.

    """

    # Pull all other filepaths that need editing
    for keys in config.keys():
        outer_loop = keys
        for keys in config[outer_loop].keys():
            second_loop = keys
            if hasattr(config[outer_loop][second_loop], 'keys') is False:
                if second_loop == 'slopefile' or second_loop == 'file':
                    path_to_correct = config[outer_loop][second_loop]
                    config[outer_loop][second_loop] = \
                        os.path.join(input_path, path_to_correct)
            else:
                for keys in config[outer_loop][second_loop].keys():
                    third_loop = keys
                    if hasattr(config[outer_loop][second_loop][third_loop],
                               'keys') is False:
                        if third_loop == 'file' or third_loop == 'filepath':
                            path_to_correct = \
                                config[outer_loop][second_loop][third_loop]
                            config[outer_loop][second_loop][third_loop] = \
                                os.path.join(input_path, path_to_correct)
                    else:
                        for keys in config[outer_loop][second_loop][third_loop].keys():
                            fourth_loop = keys
                            if hasattr(config[outer_loop][second_loop][third_loop][fourth_loop], 'keys') is False:
                                if fourth_loop == 'file' or fourth_loop == 'filepath':
                                    path_to_correct = config[outer_loop][second_loop][third_loop][fourth_loop]
                                    config[outer_loop][second_loop][third_loop][fourth_loop] = os.path.join(
                                        input_path, path_to_correct)
                            else:
                                for keys in config[outer_loop][second_loop][third_loop][fourth_loop].keys():
                                    fifth_loop = keys
                                    if fifth_loop == 'file' or fifth_loop == 'filepath':
                                        path_to_correct = config[outer_loop][second_loop][third_loop][fourth_loop][fifth_loop]
                                        config[outer_loop][second_loop][third_loop][fourth_loop][fifth_loop] = os.path.join(
                                            input_path, path_to_correct)

    return config


def validate(configfile, inputfilepath=None):
    '''
    Return a validated config object.

    Args:
        configfile (str): Config file to validate.
        inputfilepath (str): Path to input file.

    Returns:
        A validated ConfigObj object or a dictionary of which
        section/parameters failed validation.
    '''
    thispath = os.path.dirname(os.path.abspath(__file__))
    configspec = os.path.join(thispath, 'configspec.ini')
    config = ConfigObj(configfile, configspec=configspec)
    if inputfilepath is not None:
        config = correct_config_filepaths(config)
    validator = __getCustomValidator()
    result = config.validate(validator, preserve_errors=True)
    if result is True:
        return config
    else:
        errormsg = __filterResults(result)
        raise VdtTypeError(errormsg)

    return config
