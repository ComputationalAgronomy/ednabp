import tempfile
from unittest.mock import Mock

import pytest

from ednabp.common.config import Config


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
