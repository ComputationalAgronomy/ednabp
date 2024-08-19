import logging

# FORMAT = "%(asctime)s [%(levelname)s: %(message)s] (%(filename)s:%(lineno)d)"
# TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
FORMAT = "==LOG== %(message)s (%(filename)s:%(lineno)d)"
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter(FORMAT) #, TIME_FORMAT)

sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
logger.addHandler(sh)

def _get_file_handler(log_path):
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    return fh
