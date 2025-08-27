from unittest.mock import Mock

import pytest

from ednabp.bp.step_build.function_runner import FunctionRunner
from ednabp.common.config import Config


class TestFunctionRunner:
    @pytest.fixture
    def mock_config(self):
        config = Mock(spec=Config)
        config.verbose = True
        config.dry = False
        config.logger = Mock()
        return config

    @pytest.fixture
    def dry_config(self):
        config = Mock(spec=Config)
        config.verbose = True
        config.dry = True
        config.logger = Mock()
        return config

    def test_runner_init(self, mock_config):
        def test_function():
            return "test_result"

        runner = FunctionRunner("test_function", test_function, mock_config)

        assert runner.prog_name == "test_function"
        assert runner.function == test_function
        assert runner.config == mock_config
        assert runner.message == "Program: test_function."
        assert runner.capture_output is None

    def test_run_success(self, mock_config):
        expected_result = "success_result"

        def test_function():
            return expected_result

        runner = FunctionRunner("test_function", test_function, mock_config)
        result = runner.run()

        assert result is True
        assert runner.capture_output == expected_result
        mock_config.logger.info.assert_any_call("Program: test_function.")
        mock_config.logger.info.assert_any_call("COMPLETE: test_function.")

    def test_run_with_side_effects(self, mock_config):
        side_effect_tracker = []

        def test_function():
            side_effect_tracker.append("executed")
            return "done"

        runner = FunctionRunner(
            "side_effect_function", test_function, mock_config
        )
        result = runner.run()

        assert result is True
        assert side_effect_tracker == ["executed"]
        assert runner.capture_output == "done"

    def test_run_dry(self, dry_config):
        execution_tracker = []

        def test_function():
            execution_tracker.append("executed")
            return "result"

        runner = FunctionRunner("dry_function", test_function, dry_config)
        result = runner.run()

        assert result is None
        assert execution_tracker == []
        dry_config.logger.info.assert_called_with("Program: dry_function.")

    def test_run_function_value_error(self, mock_config):
        def value_error_function():
            raise ValueError("Test error")

        runner = FunctionRunner(
            "value_error_function", value_error_function, mock_config
        )
        result = runner.run()

        assert result is False
        mock_config.logger.error.assert_called()
        error_call = mock_config.logger.error.call_args[0][0]
        assert "FAIL: value_error_function" in error_call
        assert "Exception: Test error" in error_call

    def test_run_function_runtime_error(self, mock_config):
        def runtime_error_function():
            raise RuntimeError("Runtime error occurred")

        runner = FunctionRunner(
            "runtime_error_function", runtime_error_function, mock_config
        )
        result = runner.run()

        assert result is False
        mock_config.logger.error.assert_called()
        error_call = mock_config.logger.error.call_args[0][0]
        assert "FAIL: runtime_error_function" in error_call
        assert "Exception: Runtime error occurred" in error_call

    def test_run_function_type_error(self, mock_config):
        def type_error_function():
            return "foobar" + 42

        runner = FunctionRunner(
            "type_error_function", type_error_function, mock_config
        )
        result = runner.run()

        assert result is False
        mock_config.logger.error.assert_called()
        error_call = mock_config.logger.error.call_args[0][0]
        assert "FAIL: type_error_function" in error_call
        assert (
            """Exception: can only concatenate str (not "int") to str"""
            in error_call
        )
