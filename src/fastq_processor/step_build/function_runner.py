from fastq_processor.step_build import runner
from fastq_processor.step_build.stage_config import StageConfig


class FunctionRunner(runner.Runner):
    def __init__(self, prog_name, function, config):
        super().__init__(prog_name, config)
        self.function = function

    def run(self) -> bool:
        """
        Executes the function with the provided arguments.

        :returns: True if the execution is successful, False otherwise.
        """
        if self.verbose:
            self.logger.info(self.message)
        
        if not self.dry:
            try:
                self.function()
                self.logger.info(f"COMPLETE: {self.prog_name}.")
                return True
            except Exception as e:
                self.logger.error(f"FAIL: {self.prog_name}. Exception: {e.stderr}")
                return False