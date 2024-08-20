import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import umap

from analysis_toolkit.runner_build import base_runner


class HeatmapRunner(base_runner.AbundanceRunner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)

    def run_write(self,
            write_type: str = "abundance",
            taxa_level: str = "species",
            save_dir: str = ".",
            normalize: bool = False,
            sample_id_list: list[str] = []
        ):
        return super().run_write(
            write_type=write_type,
            taxa_level=taxa_level,
            save_dir=save_dir,
            normalize=normalize,
            sample_id_list=sample_id_list
        )

    @base_runner.log_execution("Plot heatmap", "plot_heatmap.log")
    def run_plot(self,
            csv_path: str,
            x_categories: list[str],
            save_dir: str = None,
            dereplicate: bool = False,
            colorscale: str = "tempo"
        ):
        """
        Plot a heatmap of the data.

        :param csv_path: The path to the CSV file containing the data.
        :param x_categories: The categories to use for the x-axis. Can be a single string or a list of strings.
        :param save_dir: The directory to save the .HTML file. If not provided, the output will not be saved. Default is None.
        :param dereplicate: Whether to dereplicate by "Sample". Default is False.
        """
        self.df = pd.read_csv(csv_path)

        if dereplicate:
            self.df = self.df.drop("Sample", axis=1)

        self.pivot_table = pd.pivot_table(self.df,
            values="Counts",
            index=self.df.columns[0],
            columns=["Site", "Year", "Month", "Sample"],
            aggfunc="mean"
        )
        self.pivot_table = self.pivot_table.fillna(0)

        self._sort_index()

        self._load_heatmap_x(x_categories)

        self.fig = go.Figure(
            data=go.Heatmap(
                z=np.array(self.pivot_table)[self.s_index,:],
                y=self.pivot_table.index[self.s_index],
                x=self.x,
                colorscale=colorscale
            )
        )
        self.fig.update_layout(yaxis_title=None, width=1000, height=600)
        self.fig.show()

        if save_dir:
            self._save_html("Heatmap", save_dir, os.path.basename(csv_path).split(".")[0])

        self.analysis_type = "Plot heatmap"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "csv_path": csv_path,
                "x_categories": x_categories,
                "dereplicate": dereplicate,
            }
        )

    def _load_heatmap_x(self, x_categories):
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

    def _sort_index(self):
        self.s_index = umap.UMAP(n_components=1, n_neighbors=15).fit(np.array(self.pivot_table)).embedding_
        self.s_index = np.argsort(self.s_index[:,0])