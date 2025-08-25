import logging
from unittest.mock import Mock, patch

import pytest

from ednabp.common.config import Config
from ednabp.common.default_settings import SETTINGS


class TestConfig:
    """Test suite for Config class functionality."""

    @pytest.fixture
    def mock_logger(self):
        return Mock(spec=logging.Logger)

    @pytest.fixture
    def config(self, mock_logger):
        return Config(logger=mock_logger, verbose=False, dry=False)

    def test_config_init_defaults(self):
        config = Config()

        assert config.verbose == SETTINGS["config_basic"]["verbose"]
        assert config.dry == SETTINGS["config_basic"]["dry"]
        assert config.logger == SETTINGS["config_basic"]["logger"]
        assert config.config_categories == ["basic"]

    def test_config_initialization_custom_values(self, mock_logger):
        config = Config(logger=mock_logger, verbose=True, dry=True)

        assert config.verbose is True
        assert config.dry is True
        assert config.logger == mock_logger
        assert config.config_categories == ["basic"]

    def test_verbose_setter_enables_verbose(self, config):
        config.verbose = True

        config.logger.setLevel.assert_called_with("INFO")

    def test_verbose_setter_disables_verbose(self, config):
        config.verbose = False

        config.logger.setLevel.assert_called_with("WARNING")

    def test_get_basic_config(self, config):
        """Test get_basic_config returns correct dictionary."""
        result = config.get_basic_config()

        expected = {
            "logger": config.logger,
            "verbose": config.verbose,
            "dry": config.dry,
        }
        assert result == expected

    def test_add_machine_config(self, config):
        config.add_machine_config(n_cpu=4, memory=16)

        assert "machine" in config.config_categories
        assert config.n_cpu == 4
        assert config.memory == 16

    def test_get_machine_config_added(self, config):
        config.add_machine_config(n_cpu=8, memory=32)

        result = config.get_machine_config()

        expected = {"n_cpu": 8, "memory": 32}
        assert result == expected

    def test_get_machine_config_not_added(self, config):
        result = config.get_machine_config()

        config.logger.warning.assert_called_once()
        assert result is None

    def test_add_iqtree_config(self, config):
        config.add_iqtree_config(model="GTR+G", boostrap=1000, overwrite=True)

        assert "iqtree" in config.config_categories
        assert config.iqtree_model == "GTR+G"
        assert config.iqtree_boostrap == 1000
        assert config.iqtree_overwrite is True

    def test_get_iqtree_config_added(self, config):
        config.add_machine_config(n_cpu=4, memory=8)
        config.add_iqtree_config(model="GTR+G", boostrap=1000, overwrite=True)

        result = config.get_iqtree_config()

        expected = {
            "threads": 4,
            "model": "GTR+G",
            "boostrap": 1000,
            "overwrite": True,
        }
        assert result == expected

    def test_get_iqtree_config_not_added_iqtree(self, config):
        config.add_machine_config(n_cpu=4, memory=8)
        result = config.get_iqtree_config()

        config.logger.warning.assert_called_once()
        assert result is None

    def test_get_iqtree_config_not_added_machine(self, config):
        config.add_iqtree_config(model="GTR+G", boostrap=1000, overwrite=True)
        result = config.get_iqtree_config()

        config.logger.warning.assert_called_once()
        assert result is None

    def test_add_cluster_config(self, config):
        reducer_kwargs = {"n_components": 2}
        clusterer_kwargs = {"min_cluster_size": 5}

        config.add_cluster_config(
            reducer_kwargs, clusterer_kwargs, encode="onehot"
        )

        assert "cluster" in config.config_categories
        assert config.cluster_reducer_kwargs == reducer_kwargs
        assert config.cluster_clusterer_kwargs == clusterer_kwargs
        assert config.cluster_encode == "onehot"

    def test_get_cluster_config_added(self, config):
        reducer_kwargs = {"n_components": 2}
        clusterer_kwargs = {"min_cluster_size": 5}

        config.add_cluster_config(
            reducer_kwargs, clusterer_kwargs, encode="onehot"
        )

        result = config.get_cluster_config()

        expected = {
            "reducer_kwargs": reducer_kwargs,
            "clusterer_kwargs": clusterer_kwargs,
            "encode": "onehot",
        }
        assert result == expected

    def test_get_cluster_config_not_added(self, config):
        result = config.get_cluster_config()

        config.logger.warning.assert_called_once()
        assert result is None

    def test_add_plot_config(self, config):
        config.add_plot_config(
            show_plot=True, save_dir="/tmp", overwrite=False
        )

        assert "plot" in config.config_categories
        assert config.plot_show_plot is True
        assert config.plot_save_dir == "/tmp"
        assert config.plot_overwrite is False

    def test_get_plot_config_added(self, config):
        config.add_plot_config(
            show_plot=True, save_dir="/tmp", overwrite=False
        )

        result = config.get_plot_config()

        expected = {
            "show_plot": True,
            "save_dir": "/tmp",
            "overwrite": False,
        }
        assert result == expected

    def test_get_plot_config_not_added(self, config):
        result = config.get_plot_config()

        config.logger.warning.assert_called_once()
        assert result is None

    def test_setattr_with_uninitialized_logger(self):
        config = Config.__new__(Config)

        # This should not raise an AttributeError
        config.verbose = True

        assert config.verbose is True
