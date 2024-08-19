from abc import ABC, abstractmethod

from fastq_processor.step_build import stage_config


class Runner(ABC):

    """
    Abstract base class for defining a runner that executes a program.

    Attributes:
        prog_name (str): The name of the program to run.
        message (str): Log message indicating the program being run.
        capture_output (subprocess.CompletedProcess): Captures the output of the subprocess.
        verbose (bool): Indicates if verbose logging is enabled.
        dry (bool): Indicates if the runner should perform a dry run (no actual execution).
        logger (Logger): Logger object for logging messages.
    """

    MSG_LOG = "==LOG=="
    MSG_DEBUG = "===DEBUG==="

    def __init__(self, prog_name: str, config: stage_config.StageConfig):
        self.prog_name = prog_name  # for debug/naming only
        self.message = f"Program: {self.prog_name}."
        self.capture_output = None
        try:
            self.verbose = config.verbose
            self.dry = config.dry
            self.logger = config.logger
        except AttributeError as e:
            print(e)


    @abstractmethod
    def run(self):
        """
        Return True/False for success/failure.
        """
        pass
