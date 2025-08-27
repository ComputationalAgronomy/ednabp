from unittest.mock import Mock

import pytest

from ednabp.common.config import Config


@pytest.fixture(scope="session")
def mock_config():
    config = Mock(spec=Config)
    config.verbose = True
    config.dry = False
    config.logger = Mock()
    return config


@pytest.fixture(scope="session")
def dry_config():
    config = Mock(spec=Config)
    config.verbose = True
    config.dry = True
    config.logger = Mock()
    return config
