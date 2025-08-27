import os
import subprocess
import tempfile
from unittest.mock import Mock, patch

import pytest

from ednabp.bp.step_build.subproces_runner import (
    RedirectOutputRunner,
    SubprocessRunner,
)
from ednabp.common.config import Config


class TestSubprocessRunner:
    @pytest.fixture
    def subprocess_runner(self, mock_config):
        return SubprocessRunner("test_subprocess", "echo test", mock_config)

    def test_runner_init_default(self, mock_config):
        runner = SubprocessRunner("test_subprocess", "echo test", mock_config)

        assert runner.prog_name == "test_subprocess"
        assert runner.command == "echo test"
        assert runner.config == mock_config
        assert runner.shell is False
        assert runner.message == "Program: test_subprocess."
        assert runner.capture_output is None

    def test_runner_init_custom(self, mock_config):
        runner = SubprocessRunner(
            "test_subprocess", "echo test", mock_config, shell=True
        )
        assert runner.shell is True

    @patch("subprocess.run")
    def test_run_success(self, mock_subprocess_run, subprocess_runner):
        mock_result = Mock()
        mock_result.returncode = 0
        mock_subprocess_run.return_value = mock_result

        result = subprocess_runner.run()

        assert result is True
        subprocess_runner.config.logger.info.assert_any_call(
            "Program: test_subprocess."
        )
        subprocess_runner.config.logger.info.assert_any_call(
            "COMPLETE: test_subprocess."
        )
        mock_subprocess_run.assert_called_once()

    @patch("subprocess.run")
    def test_run_shell_false(self, mock_subprocess_run, mock_config):
        mock_result = Mock()
        mock_subprocess_run.return_value = mock_result

        runner = SubprocessRunner(
            "test_subprocess", "echo test", mock_config, shell=False
        )
        runner.run()

        call_args = mock_subprocess_run.call_args[0][0]
        assert isinstance(call_args, list)
        assert call_args == ["echo", "test"]

    @patch("subprocess.run")
    def test_run_shell_true(self, mock_subprocess_run, mock_config):
        mock_result = Mock()
        mock_subprocess_run.return_value = mock_result

        runner = SubprocessRunner(
            "test_subprocess", "echo test", mock_config, shell=True
        )
        runner.run()

        call_args = mock_subprocess_run.call_args[0][0]
        assert call_args == "echo test"

        call_kwargs = mock_subprocess_run.call_args[1]
        assert call_kwargs["shell"] is True

    @patch("subprocess.run")
    def test_run_called_process_error(
        self, mock_subprocess_run, subprocess_runner
    ):
        error = subprocess.CalledProcessError(
            1, "echo", stderr="Command failed"
        )
        mock_subprocess_run.side_effect = error

        result = subprocess_runner.run()

        assert result is False
        subprocess_runner.config.logger.error.assert_called_with(
            "FAIL: test_subprocess. SubprocessError: Command failed."
        )

    @patch("subprocess.run")
    def test_run_file_not_found_error(
        self, mock_subprocess_run, subprocess_runner
    ):
        error = FileNotFoundError("Command not found")
        error.strerror = "No such file or directory"
        mock_subprocess_run.side_effect = error

        result = subprocess_runner.run()

        assert result is False
        subprocess_runner.config.logger.error.assert_called_with(
            "FAIL: test_subprocess. FileNotFoundError: No such file or directory."
        )

    @patch("subprocess.run")
    def test_run_other_exception(self, mock_subprocess_run, subprocess_runner):
        error = RuntimeError("Unexpected error")
        mock_subprocess_run.side_effect = error

        result = subprocess_runner.run()

        assert result is False
        subprocess_runner.config.logger.error.assert_called_with(
            "FAIL: test_subprocess. Other Exception: Unexpected error."
        )

    def test_run_dry(self, dry_config):
        runner = SubprocessRunner("test_subprocess", "echo test", dry_config)

        result = runner.run()

        assert result is None
        dry_config.logger.info.assert_called_with("Program: test_subprocess.")

    @patch("subprocess.run")
    def test_capture_output_stored(
        self, mock_subprocess_run, subprocess_runner
    ):
        mock_result = Mock()
        mock_result.stdout = "test output"
        mock_result.stderr = "test error"
        mock_subprocess_run.return_value = mock_result

        subprocess_runner.run()

        assert subprocess_runner.capture_output == mock_result
        assert subprocess_runner.capture_output.stdout == "test output"
        assert subprocess_runner.capture_output.stderr == "test error"

    @patch("shlex.split")
    @patch("subprocess.run")
    def test_command_splitting_windows(
        self, mock_subprocess_run, mock_shlex_split, mock_config
    ):
        mock_result = Mock()
        mock_subprocess_run.return_value = mock_result
        mock_shlex_split.return_value = ["echo", "test"]

        with patch("sys.platform", "win32"):
            runner = SubprocessRunner(
                "test_subprocess", "echo test", mock_config, shell=False
            )
            runner.run()

        mock_shlex_split.assert_called_with("echo test", posix=False)

    @patch("shlex.split")
    @patch("subprocess.run")
    def test_command_splitting_linux(
        self, mock_subprocess_run, mock_shlex_split, mock_config
    ):
        mock_result = Mock()
        mock_subprocess_run.return_value = mock_result
        mock_shlex_split.return_value = ["echo", "test"]

        with patch("sys.platform", "linux"):
            runner = SubprocessRunner(
                "test_subprocess", "echo test", mock_config, shell=False
            )
            runner.run()

        mock_shlex_split.assert_called_with("echo test", posix=True)

    @patch("subprocess.run")
    def test_quote_removal(self, mock_subprocess_run, mock_config):
        mock_result = Mock()
        mock_subprocess_run.return_value = mock_result

        runner = SubprocessRunner(
            "test_subprocess", 'echo "test world"', mock_config, shell=False
        )

        with patch("shlex.split", return_value=["echo", '"test world"']):
            runner.run()

        call_args = mock_subprocess_run.call_args[0][0]
        assert call_args == ["echo", "test world"]


class TestRedirectOutputRunner:
    @pytest.fixture
    def mock_subprocess_runner(self, mock_config):
        runner = Mock(spec=SubprocessRunner)
        runner.capture_output = Mock()
        runner.capture_output.stdout = "test stdout"
        runner.capture_output.stderr = "test stderr"
        return runner

    @pytest.fixture
    def tmp_files(self):
        with (
            tempfile.NamedTemporaryFile(mode="w", delete=False) as stdout_file,
            tempfile.NamedTemporaryFile(mode="w", delete=False) as stderr_file,
        ):
            stdout_path = stdout_file.name
            stderr_path = stderr_file.name

        yield stdout_path, stderr_path

        try:
            os.unlink(stdout_path)
            os.unlink(stderr_path)
        except FileNotFoundError:
            pass

    def test_runner_init(self, mock_config, mock_subprocess_runner, tmp_files):
        stdout_file, stderr_file = tmp_files

        runner = RedirectOutputRunner(
            "test_redirect",
            mock_subprocess_runner,
            stdout_file,
            stderr_file,
            mock_config,
        )

        assert runner.prog_name == "test_redirect"
        assert runner.runner == mock_subprocess_runner
        assert runner.stdout_file == stdout_file
        assert runner.stderr_file == stderr_file
        assert runner.config == mock_config
        assert runner.message == "RedirectOutput: test_redirect."

    def test_invalid_subprocess_runner(self, mock_config, tmp_files):
        stdout_file, stderr_file = tmp_files
        invalid_runner = "invalid_runner"

        with pytest.raises(TypeError, match="Invalid instance type"):
            RedirectOutputRunner(
                "test", invalid_runner, stdout_file, stderr_file, mock_config
            )

    def test_run_success(self, mock_config, mock_subprocess_runner, tmp_files):
        stdout_file, stderr_file = tmp_files

        runner = RedirectOutputRunner(
            "test_redirect",
            mock_subprocess_runner,
            stdout_file,
            stderr_file,
            mock_config,
        )

        result = runner.run()

        assert result is True
        mock_config.logger.info.assert_any_call(
            "RedirectOutput: test_redirect."
        )
        mock_config.logger.info.assert_any_call("COMPLETE: test_redirect.")

        with open(stdout_file) as out_f, open(stderr_file) as err_f:
            assert out_f.read() == "test stdout"
            assert err_f.read() == "test stderr"

    def test_run_no_stdout(
        self, mock_config, mock_subprocess_runner, tmp_files
    ):
        stdout_file, stderr_file = tmp_files
        mock_subprocess_runner.capture_output.stdout = None

        runner = RedirectOutputRunner(
            "test_redirect",
            mock_subprocess_runner,
            stdout_file,
            stderr_file,
            mock_config,
        )

        result = runner.run()

        assert result is True

        with open(stdout_file) as out_f, open(stderr_file) as err_f:
            assert out_f.read() == ""
            assert err_f.read() == "test stderr"

    def test_run_no_stderr(
        self, mock_config, mock_subprocess_runner, tmp_files
    ):
        stdout_file, stderr_file = tmp_files
        mock_subprocess_runner.capture_output.stderr = None

        runner = RedirectOutputRunner(
            "test_redirect",
            mock_subprocess_runner,
            stdout_file,
            stderr_file,
            mock_config,
        )

        result = runner.run()

        assert result is True

        with open(stdout_file) as out_f, open(stderr_file) as err_f:
            assert out_f.read() == "test stdout"
            assert err_f.read() == ""

    def test_run_no_output(self, mock_config, tmp_files):
        stdout_file, stderr_file = tmp_files
        mock_subprocess_runner = Mock(spec=SubprocessRunner)
        mock_subprocess_runner.capture_output = None

        runner = RedirectOutputRunner(
            "test_redirect",
            mock_subprocess_runner,
            stdout_file,
            stderr_file,
            mock_config,
        )

        result = runner.run()

        assert result is False
        mock_config.logger.error.assert_called()
        error_call = mock_config.logger.error.call_args[0][0]
        assert "FAIL: test_redirect. Error:" in error_call

    def test_run_dry(self, mock_subprocess_runner, tmp_files):
        stdout_file, stderr_file = tmp_files
        dry_config = Mock(spec=Config)
        dry_config.verbose = True
        dry_config.dry = True
        dry_config.logger = Mock()

        runner = RedirectOutputRunner(
            "test_redirect",
            mock_subprocess_runner,
            stdout_file,
            stderr_file,
            dry_config,
        )

        result = runner.run()

        assert result is None
        dry_config.logger.info.assert_called_with(
            "RedirectOutput: test_redirect."
        )

        with open(stdout_file) as f:
            assert f.read() == ""

    def test_file_creation(self, mock_config, mock_subprocess_runner):
        with tempfile.TemporaryDirectory() as temp_dir:
            stdout_file = os.path.join(temp_dir, "stdout.txt")
            stderr_file = os.path.join(temp_dir, "stderr.txt")

            runner = RedirectOutputRunner(
                "test_redirect",
                mock_subprocess_runner,
                stdout_file,
                stderr_file,
                mock_config,
            )

            result = runner.run()

            assert result is True
            assert os.path.exists(stdout_file)
            assert os.path.exists(stderr_file)
