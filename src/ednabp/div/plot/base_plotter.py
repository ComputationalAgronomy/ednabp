import os
from abc import ABC, abstractmethod

import pandas as pd

from ...common import config

VALUE_COLUMN = "values"


class Plotter(ABC):
    @abstractmethod
    def __init__(self, verbose, show_plot, save_dir, overwrite):
        self.config = config.Config(verbose=verbose)
        self.config.add_plot_config(show_plot, save_dir, overwrite)

    def read_df(self, df):
        if isinstance(df, pd.DataFrame):
            self.df = df
        elif isinstance(df, str):
            self.df = pd.read_csv(df)
        else:
            raise ValueError(
                "df must be a pandas DataFrame or a path to a csv file"
            )

    @abstractmethod
    def plot(self):
        pass

    def show_and_save(
        self,
        fig,
        prefix,
    ):
        config = self.config.get_plot_config()
        save_dir = config["save_dir"]

        if config["show_plot"]:
            fig.show()

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            html_path = os.path.join(save_dir, f"{prefix}.html")
            if not os.path.exists(html_path) or config["overwrite"]:
                fig.write_html(html_path)
                self.config.logger.info(f"Saved plot to {html_path}")
