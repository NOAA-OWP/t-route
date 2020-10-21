import logging

def set_logger(LOG,verbose,debuglevel,log_writer,log_file):
    if verbose:
        LOG.setLevel(logging.INFO)
    elif debuglevel == 1:
        LOG.setLevel(logging.DEBUG)
    elif debuglevel == 2:
        LOG.setLevel(logging.WARNING)
    elif debuglevel == 3:
        LOG.setLevel(logging.CRITICAL)
    else:
        LOG.setLevel(logging.NOTSET)
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        filename=log_file,
        filemode=log_writer,
    )
    return LOG

    