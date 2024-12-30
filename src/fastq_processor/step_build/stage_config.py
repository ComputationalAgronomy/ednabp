import sys


class StageConfig:
    """
    Configuration class for stage settings.

    Attributes:
        verbose (bool): Indicates if verbose logging is enabled.
        dry (bool): Indicates if the stage should perform a dry run.
        logger (Logger): Logger object for logging messages.
        n_cpu (int): Number of CPU cores allocated for the stage.
        memory (int): Amount of memory (in GB) allocated for the stage.
    """
    def __init__(self,
                 verbose: bool = False,
                 dry: bool = False,
                 logger = sys.stdout,
                 n_cpu: int = 1,
                 memory: int = 8
                 ):
        self.verbose = verbose
        self.dry = dry
        self.logger = logger
        self.n_cpu = n_cpu
        self.memory = memory

    def get_machine_info(self) -> dict[str, int]:
        """
        Retrieves the machine information including CPU cores and memory.

        :returns: dict: A dictionary containing the number of CPU cores and amount of memory.
        """
        return {"n_cpu": self.n_cpu,
                "memory": self.memory
                }

    def get_basic_configuration(self) -> dict:
        """
        Retrieves the basic configuration settings.

        :returns: dict: A dictionary containing the verbose, dry run, and logger settings.
        """
        return {"verbose": self.verbose,
                "dry_run": self.dry,
                "logger": self.logger
                }

    # TODO(SW): Use these to help you organise all these parameters. dict[str:dict]
    def get_usearch_configuration(self) -> dict:
        return {}

    def get_denoise_configuration(self) -> dict:
        return {}