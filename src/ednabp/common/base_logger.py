import logging
import os
import sys

# FORMAT = "%(asctime)s [%(levelname)s: %(message)s] (%(filename)s:%(lineno)d)"
# TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
FORMAT = "==%(levelname)s== %(message)s"  # (%(filename)s:%(lineno)d)"
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter(FORMAT)  # , TIME_FORMAT)

sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
logger.addHandler(sh)


def get_file_handler(log_path):
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    return fh


def close_file_handler(logger):
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            logger.removeHandler(handler)
            handler.close()


def prog_log(prog_name: str, log_file: str | None = None):
    def decorator(func):
        def wrapper(*args, **kwargs):
            if (
                log_file is not None
                and "save_dir" in kwargs
                and os.path.exists(kwargs["save_dir"])
            ):
                fh = get_file_handler(
                    os.path.join(kwargs["save_dir"], log_file)
                )
                logger.addHandler(fh)
                write_log = True
            else:
                write_log = False
            try:
                logger.info(f"Program: {prog_name}")
                result = func(*args, **kwargs)
                logger.info(f"COMPLETE: {prog_name}")
                return result
            except Exception as e:
                logger.error(f"FAIL: {prog_name}. {type(e).__name__}: {e}")
                sys.exit(0)
            finally:
                if write_log:
                    logger.removeHandler(fh)

        return wrapper

    return decorator
