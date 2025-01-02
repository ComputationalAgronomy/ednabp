from abc import ABC, abstractmethod
import os

from . import base_logger

class DMPlotter(ABC):
    SAMPLE_ID_COLUMN = "Sample_id"
    logger = base_logger.logger

    def __init__(self):
        pass

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

    @base_logger.prog_log("Display and save plot (if 'save_dir' provided)")
    def _display_and_save_plot(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            self._save_plot(save_dir, overwrite)

    @abstractmethod
    def _save_plot(self):
        pass