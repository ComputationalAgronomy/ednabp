import os
import tempfile
from unittest.mock import patch

import pytest

from ednabp.bp.step_build.function_runner import FunctionRunner
from ednabp.bp.step_build.stage_builder import StageBuilder
from ednabp.bp.step_build.subproces_runner import SubprocessRunner


class ImplementedStageBuilder(StageBuilder):
    def __init__(self, config, in_dir, out_dir):
        super().__init__("test_stage", config, in_dir, out_dir)
        self.in_suffix = "_test.txt"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")


class TestStageBuilder:
    @pytest.fixture
    def tmp_dirs(self):
        with (
            tempfile.TemporaryDirectory() as in_dir,
            tempfile.TemporaryDirectory() as out_dir,
        ):
            yield in_dir, out_dir

    @pytest.fixture
    def stage_builder(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return ImplementedStageBuilder(mock_config, in_dir, out_dir)

    def test_stage_builder_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = ImplementedStageBuilder(mock_config, in_dir, out_dir)

        assert stage.heading == "test_stage"
        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.runners == []
        assert stage.output == []
        assert os.path.exists(out_dir)

    def test_add_stage(self, stage_builder):
        stage_builder.add_stage("test_stage", "echo test", shell=False)

        assert len(stage_builder.runners) == 1
        runner = stage_builder.runners[0]
        assert isinstance(runner, SubprocessRunner)
        assert runner.prog_name == "test_stage"
        assert runner.command == "echo test"
        assert runner.shell is False

    def test_add_stage_shell_true(self, stage_builder):
        stage_builder.add_stage("test_stage", "echo test", shell=True)

        runner = stage_builder.runners[0]
        assert runner.shell is True

    @patch("ednabp.bp.step_build.subproces_runner.RedirectOutputRunner")
    def test_add_stage_output_to_file(
        self, mock_redirect_runner, stage_builder
    ):
        stage_builder.add_stage("test_stage", "echo test")

        stage_builder.add_stage_output_to_file(
            "redirect_test", 0, "output.txt", "error.txt"
        )

        assert len(stage_builder.runners) == 2
        mock_redirect_runner.assert_called_once_with(
            "redirect_test",
            stage_builder.runners[0],
            "output.txt",
            "error.txt",
            stage_builder.config,
        )

    def test_add_stage_function(self, stage_builder):
        def test_function():
            return "test_result"

        stage_builder.add_stage_function("test_function", test_function)

        assert len(stage_builder.runners) == 1
        runner = stage_builder.runners[0]
        assert isinstance(runner, FunctionRunner)
        assert runner.prog_name == "test_function"
        assert runner.function == test_function

    def test_setup_infile_creation(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = ImplementedStageBuilder(mock_config, in_dir, out_dir)

        stage.setup("test_infile")

        expected_infile = os.path.join(in_dir, "test_infile_test.txt")
        assert stage.infile == expected_infile

    def test_check_infile(self, stage_builder, tmp_dirs):
        in_dir, _ = tmp_dirs
        test_file = os.path.join(in_dir, "test_test.txt")

        with open(test_file, "w") as f:
            f.write("test content")

        stage_builder.setup("test")
        stage_builder.check_infile()

    def test_check_infile_file_not_found_error(self, stage_builder):
        stage_builder.setup("nonexistent")

        with pytest.raises(
            FileNotFoundError, match="nonexistent_test.txt not found"
        ):
            stage_builder.check_infile()

    def test_check_infile_attribute_error(self, stage_builder):
        with pytest.raises(AttributeError):
            stage_builder.check_infile()

    def test_check_outdir(self, mock_config):
        with tempfile.TemporaryDirectory() as temp_dir:
            in_dir = temp_dir
            out_dir = os.path.join(temp_dir, "new_output_dir")
            assert not os.path.exists(out_dir)

            _stage = ImplementedStageBuilder(mock_config, in_dir, out_dir)
            assert os.path.exists(out_dir)

    def test_summary(self, stage_builder):
        stage_builder.add_stage("program1", "echo test")
        stage_builder.add_stage_function("function1", lambda: None)

        summary = stage_builder.summary()

        expected = [
            "Step 0: Program: program1.",
            "Step 1: Program: function1.",
        ]
        assert summary == expected

    @patch("ednabp.bp.step_build.subproces_runner.SubprocessRunner.run")
    @patch("ednabp.bp.step_build.function_runner.FunctionRunner.run")
    def test_run_success(
        self, mock_func_run, mock_subprocess_run, stage_builder
    ):
        mock_subprocess_run.return_value = True
        mock_func_run.return_value = True

        stage_builder.add_stage("program1", "echo test")
        stage_builder.add_stage_function("function1", lambda: None)

        result = stage_builder.run()

        assert result is True
        assert stage_builder.output == [True, True]
        assert len(stage_builder.runners) == 0
        stage_builder.config.logger.info.assert_called_with(
            "Running: test_stage"
        )

    @patch("ednabp.bp.step_build.subproces_runner.SubprocessRunner.run")
    def test_run_failure(self, mock_run, stage_builder):
        mock_run.side_effect = [True, False]

        stage_builder.add_stage("program1", "echo test")
        stage_builder.add_stage("program2", "false")

        result = stage_builder.run()

        assert result is False
        assert stage_builder.output == [True, False]

    @patch("os.makedirs")
    def test_check_outdir_existed(self, mock_makedirs, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        ImplementedStageBuilder(mock_config, in_dir, out_dir)

        mock_makedirs.assert_called_with(out_dir, exist_ok=True)
