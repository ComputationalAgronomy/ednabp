import shlex
import subprocess
import sys

from fastq_processor.step_build import runner


class SubprocessRunner(runner.Runner):

    """
    Runner class for executing a command as a subprocess.

    Attributes:
        command (str): The command to execute.
        shell (bool): Whether to execute the command through the shell.
    """

    def __init__(self, prog_name: str, command: str, config, shell: bool = False):
        super().__init__(prog_name, config)

        self.command = command
        self.shell = shell

    def run(self) -> bool:
        """
        Executes the command as a subprocess.

        :returns: True if the execution is successful, False otherwise.
        """
        if self.verbose:
            self.logger.info(self.message)

        if self.shell:
            args = self.command
        else:
            args = shlex.split(self.command, posix="win" not in sys.platform)
            args = [arg.replace("\"", "") for arg in args]

        if not self.dry:
            try:
                self.capture_output = subprocess.run(args,
                                                   capture_output=False,
                                                   shell=self.shell,
                                                   check=True,
                                                   text=True,
                                                   stdout=subprocess.PIPE,
                                                   stderr=subprocess.PIPE)
                self.logger.info(f"COMPLETE: {self.prog_name}.")
                return True
            except subprocess.CalledProcessError as e:
                self.logger.error(f"FAIL: {self.prog_name}. SubprocessError: {e.stderr}.")
                return False
            except FileNotFoundError as e:
                self.logger.error(f"FAIL: {self.prog_name}. FileNotFoundError: {e.stderr}.")
                return False
            except Exception as e:
                self.logger.error(f"FAIL: {self.prog_name}. Other Exception: {e.stderr}.")
                return False


class RedirectOutputRunner(runner.Runner):
    """
    Runner class for redirecting the output of a subprocess to a file.

    Attributes:
        runner (SubprocessRunner): The subprocess runner whose output to redirect.
        outfile (str): The file to write the output to.
    """
    def __init__(self, prog_name: str, runner: SubprocessRunner, stdout_file: str, stderr_file: str, config):
        super().__init__(prog_name, config)

        if isinstance(runner, SubprocessRunner):
            self.runner = runner
            self.stdout_file = stdout_file
            self.stderr_file = stderr_file
            self.message = f"RedirectOutput: {self.prog_name}."
        else:
            info = type(runner)
            raise TypeError(f"Invalid instance type: {info}.")

    def run(self) -> bool:
        """
        Redirects the output of the runner to the specified file.

        :returns: True if the output redirection is successful, False otherwise.
        """
        if self.verbose:
            self.logger.info(self.message)

        if not self.dry:
            try:
                with open(self.stdout_file, "w") as out_f, open(self.stderr_file, "w") as err_f:
                    pass
                if self.runner.capture_output.stdout:
                    with open(self.stdout_file, "a") as out_f:
                        out_f.write(self.runner.capture_output.stdout)
                if self.runner.capture_output.stderr:
                    with open(self.stderr_file, "a") as err_f:
                        err_f.write(self.runner.capture_output.stderr)
                self.logger.info(f"COMPLETE: {self.prog_name}.")
                return True
            except AttributeError as e:
                self.logger.error(f"FAIL: {self.prog_name}. Error: {e}.")
                return False
