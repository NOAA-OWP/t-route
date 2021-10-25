import logging
import sys   

def log_level_set(input_parameters):
    log_level = input_parameters.get("log_level", 'DEBUG')
    logging.basicConfig(
        level=log_level,
        format='%(levelname)7s [%(filename)s:%(lineno)s - %(funcName)20s()]: %(message)s',
        stream=sys.stderr,
    )