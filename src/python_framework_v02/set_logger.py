import logging

def set_logger(LOG,verbose,debuglevel,log_writer,log_file):
    if verbose:
        LOG.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
    elif debuglevel == 1:
        LOG.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
    elif debuglevel == 2:
        LOG.setLevel(logging.WARNING)
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)
    else:
        LOG.setLevel(logging.CRITICAL)
        ch = logging.StreamHandler()
        ch.setLevel(logging.CRITICAL)
    # formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    LOG.addHandler(ch)
    # logging.basicConfig(
    #     format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    #     datefmt="%m/%d/%Y %I:%M:%S %p",
    #     filename=log_file,
    #     filemode=log_writer,
    # )
    return LOG