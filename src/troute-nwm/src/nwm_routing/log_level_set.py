import logging
import sys   

def log_level_set(input_parameters):
    '''
    Set logging level and specify logger configuration.
    
    Arguments
    ---------
    input_parameters (dict): User input logging parameters
    
    Returns
    -------
    None
    
    Notes
    -----
    In the absense of user-specified logging level, level defaults to DEGUG
    See also https://docs.python.org/3/library/logging.html
    
    '''
    log_level = input_parameters.get("log_level", 'DEBUG')
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s %(levelname)7s [%(filename)s:%(lineno)s - %(funcName)20s()]: %(message)s',
        stream=sys.stderr,
    )