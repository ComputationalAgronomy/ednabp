import tempfile
from unittest.mock import Mock

import pytest

from ednabp.common import base_logger
from ednabp.common.config import Config


@pytest.fixture(autouse=True)
def restore_logger_handlers():
    """Some pipeline tests attach mock file handlers to the shared logger;
    remove them so later tests can log normally."""
    handlers = list(base_logger.logger.handlers)
    yield
    base_logger.logger.handlers = handlers


@pytest.fixture
def mock_config():
    config = Mock(spec=Config)
    config.verbose = True
    config.dry = False
    config.logger = Mock()
    config.n_cpu = 4
    return config


@pytest.fixture
def dry_config():
    config = Mock(spec=Config)
    config.verbose = True
    config.dry = True
    config.logger = Mock()
    return config


@pytest.fixture
def tmp_dirs():
    with (
        tempfile.TemporaryDirectory() as in_dir,
        tempfile.TemporaryDirectory() as out_dir,
    ):
        yield in_dir, out_dir
