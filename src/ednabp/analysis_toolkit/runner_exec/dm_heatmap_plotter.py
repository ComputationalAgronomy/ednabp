import os
from typing import Annotated

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from ..runner_build import DMPlotter, base_logger


class HeatmapPlotter(DMPlotter):

    @base_logger.prog_log("Plot heatmap")
    def plot_heatmap(self,
            csv_path: str,
            taxa_column: str,
            metric_column: str,
            x_categories: str | Annotated[list[str], 2],
            dereplicate: bool = False,
            colorscale: str = "tempo",
            save_dir: str = None,
            overwrite: bool = False,
        ):
        """
        Plot a heatmap of the data.

        :param csv_path: The path to the CSV file containing the data.
        :param taxa_column: Column name to use for color values 
        :param metric_column: Column name to use for y-axis values
        :param x_categories: The categories to use for the x-axis. Can be a single string or a list of strings.
        :param dereplicate: Whether to dereplicate by "Sample". Default: False.
        :param colorscale: color domain within which colors are to be interpolated. Default: "tempo".
        :param save_dir: If provided, the heatmap will be saved as a .HTML file. Default: None.
        :param overwrite: Whether to overwrite existing files. Default: False.
        """
        HeatmapPlotter._load_and_validate_data(self, csv_path, taxa_column, metric_column)
        HeatmapPlotter._process_data(self, taxa_column, metric_column, dereplicate)
        HeatmapPlotter._prepare_plot_data(self, x_categories)
        HeatmapPlotter._create_plot(self, colorscale)
        HeatmapPlotter._display_and_save_plot(self, save_dir, overwrite)

    @base_logger.prog_log("Load and validate data")
    def _load_and_validate_data(self, csv_path, taxa_column, metric_column):
        try:
            self.df = pd.read_csv(csv_path)
            required_columns = {taxa_column, metric_column}
            if not required_columns.issubset(self.df.columns):
                missing = required_columns - set(self.df.columns)
                raise ValueError(f"Missing required columns: {missing}")
        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find CSV file: {csv_path}")

    @base_logger.prog_log("Create pivot table")
    def _process_data(self, taxa_column, metric_column, dereplicate):
        if dereplicate:
            self.df = self.df.drop("Sample", axis=1)
        pivot_columns = [i for i in self.df.columns.tolist() if i not in [taxa_column, metric_column]]
        self.pivot_table = pd.pivot_table(self.df,
            values=metric_column,
            index=taxa_column,
            columns=pivot_columns,
            aggfunc="mean"
        )
        self.pivot_table = self.pivot_table.fillna(0)
        self._sort_index()

    def _sort_index(self):
        import umap
        self.s_index = umap.UMAP(n_components=1, n_neighbors=15, random_state=42).fit(np.array(self.pivot_table)).embedding_
        self.s_index = np.argsort(self.s_index[:,0])

    @base_logger.prog_log("Prepare plot data")
    def _prepare_plot_data(self, x_categories):
        column_names = self.pivot_table.columns.names
        if x_categories not in column_names and not set(x_categories).issubset(set(column_names)):
            raise ValueError(f"columns must be a subset of the column names {column_names} of the pivot table")
        elif isinstance(x_categories, str):
            todrop_columns = column_names.copy()
            todrop_columns.remove(x_categories)
            self.x = self.pivot_table.columns.droplevel(todrop_columns)
        else:
            self.x = []
            for column in x_categories:
                todrop_columns = column_names.copy()
                todrop_columns.remove(column)
                self.x.append(self.pivot_table.columns.droplevel(todrop_columns))
            if len(x_categories) > 2:
                self.logger.warning("WARNING: More than two-level categorical x-axis is not yet available in Plotly yet. This is a substitute implementation that combines the first n-1 categories into the first level.")
                self.x = [["<br>".join(list(map(str, x))[::-1]) for x in zip(*self.x[:-1])], self.x[-1]]
        self.y = self.pivot_table.index[self.s_index]
        self.z = np.array(self.pivot_table)[self.s_index,:]

    @base_logger.prog_log("Create plot")
    def _create_plot(self, colorscale):
        self.fig = go.Figure(
            data=go.Heatmap(
                z=self.z,
                y=self.y,
                x=self.x,
                colorscale=colorscale
            )
        )
        self.fig.update_layout(yaxis_title=None, width=1000, height=600)

    @base_logger.prog_log("Display and save plot (if 'save_dir' provided)")
    def _display_and_save_plot(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            HeatmapPlotter._save_plot(self, save_dir, overwrite)

    def _save_plot(self, save_html_dir: str, overwrite: bool):
        fig_path = os.path.join(save_html_dir, "heatmap.html")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(f"WARNING: File already exists: {fig_path}. Stop saving.")
            return
        self.fig.write_html(fig_path)
        self.logger.info(f"Heatmap saved to: {fig_path}")