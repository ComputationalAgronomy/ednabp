from abc import ABC, abstractmethod
import os

from . import base_logger

class DMPlotter(ABC):
    SAMPLE_ID_COLUMN = "sample_id"
    logger = base_logger.logger

    def __init__(self, no_versose: bool):
        if no_versose:
            self.logger.setLevel("WARNING")

    @abstractmethod
    def _load_and_validate_data(self):
        pass

    @abstractmethod
    def _process_data(self):
        pass

    @abstractmethod
    def _prepare_plot_data(self):
        pass

    @abstractmethod
    def _create_plot(self):
        pass

    @abstractmethod
    def _display_and_save_plot(self):
        pass

    @abstractmethod
    def _save_plot(self):
        pass