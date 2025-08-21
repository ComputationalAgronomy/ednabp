from . import base_runner


class FunctionRunner(base_runner.Runner):
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
                self.logger.error(f"FAIL: {self.prog_name}. Exception: {e}")
                return False
