from abc import ABC, abstractmethod
import pandas as pd

from . import base_logger

class DMPlotter(ABC):
    SAMPLE_ID_COLUMN = "sample_id"
    logger = base_logger.logger

    def __init__(self, no_versose: bool):
        if no_versose:
            self.logger.setLevel("WARNING")
        else:
            self.logger.setLevel("INFO")

    @base_logger.prog_log("Load and validate input csv")
    def _load_and_validate_data(self, csv_path, required_columns):
        try:
            self.df = pd.read_csv(csv_path)
            if not required_columns.issubset(self.df.columns):
                missing = required_columns - set(self.df.columns)
                raise ValueError(f"Missing required columns: {missing}")
        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find CSV file: {csv_path}")

    @base_logger.prog_log("Create pivot table")
    def _create_pivot_table(self, values, index, columns, aggfunc):
        self.pivot_table = pd.pivot_table(self.df,
            values=values,
            index=index,
            columns=columns,
            aggfunc=aggfunc,
            fill_value=0
        )
        # self.pivot_table = self.pivot_table.fillna(0)

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
    def _display_and_save(self):
        pass

    @abstractmethod
    def _save_csv(self):
        pass

    @abstractmethod
    def _save_plot(self):
        pass