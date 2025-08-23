from typing import Literal, Protocol

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from . import base_plotter


# PCA(n_components=1), TSNE(n_components=1), umap.UMAP(n_components=1)
class FitTransformProtocol(Protocol):
    def fit_transform(self, X: np.ndarray) -> np.ndarray: ...


class HeatmapPlotter(base_plotter.Plotter):
    def __init__(
        self,
        df,
        values,
        index,
        columns,
        aggfunc,
        sort_method,
        verbose,
        show_plot,
        save_dir,
        overwrite,
    ):
        self.read_df(df)
        self.values = values
        self.index = index
        self.columns = columns
        self.aggfunc = aggfunc
        self.sort_method = sort_method
        self.pivot_table = None
        super().__init__(verbose, show_plot, save_dir, overwrite)

    def create_pivot_table(self):
        self.pivot_table = pd.pivot_table(
            self.df,
            values=self.values,
            index=self.index,
            columns=self.columns,
            aggfunc=self.aggfunc,
            fill_value=0,
        )

    def sort_index(self):
        data = np.array(self.pivot_table)
        if self.sort_method is None:
            self.s_index = np.arange(len(data))
        elif self.sort_method == "sum":
            self.s_index = np.argsort(data.sum(axis=1))[::-1]
        elif self.sort_method == "mean":
            self.s_index = np.argsort(data.mean(axis=1))[::-1]
        elif self.sort_method == "max":
            self.s_index = np.argsort(data.max(axis=1))[::-1]
        elif self.sort_method == "hierarchical":
            from scipy.cluster.hierarchy import dendrogram, linkage

            linkage_matrix = linkage(data, method="ward")
            dendro = dendrogram(linkage_matrix, no_plot=True)
            self.s_index = dendro["leaves"]
        elif hasattr(self.sort_method, "fit_transform"):
            embedding = self.sort_method.fit_transform(data)
            self.s_index = np.argsort(embedding[:, 0])
        else:
            self.config.logger.warning(
                f"WARNING: Unknown sort method: {self.sort_method}. Skipping the sorting procedure"
            )
            self.s_index = np.arange(len(data))

    def prepare_plot_data(self):
        column_names = self.pivot_table.columns.names
        if len(column_names) == 1:
            self.x = self.pivot_table.columns
        else:
            self.x = []
            for column in column_names:
                todrop_columns = column_names.copy()
                todrop_columns.remove(column)
                self.x.append(
                    self.pivot_table.columns.droplevel(todrop_columns)
                )
            if len(column_names) > 2:
                self.config.logger.warning(
                    "WARNING: More than two-level categorical x-axis is not yet available in Plotly yet. This is a substitute implementation that combines the first n-1 categories into the first level."
                )
                self.x = [
                    [
                        "<br>".join(list(map(str, x))[::-1])
                        for x in zip(*self.x[:-1], strict=False)
                    ],
                    self.x[-1],
                ]
        self.y = self.pivot_table.index[self.s_index]
        self.z = np.array(self.pivot_table)[self.s_index, :]

    def plot(self):
        self.create_pivot_table()
        self.sort_index()
        self.prepare_plot_data()

        fig = go.Figure(
            data=go.Heatmap(
                z=self.z,
                y=self.y,
                x=self.x,
            )
        )

        return fig


def heatmap(
    df: str | pd.DataFrame,
    values: str,
    index: str,
    columns: str | list[str],
    aggfunc: Literal["mean", "sum"] = "mean",
    sort_method: None
    | Literal["sum", "mean", "max", "hierarchical"]
    | FitTransformProtocol = None,
    verbose: bool = False,
    show_plot: bool = True,
    save_dir: str | None = None,
    overwrite: bool = False,
) -> tuple[go.Figure, HeatmapPlotter]:
    plotter = HeatmapPlotter(
        df,
        values,
        index,
        columns,
        aggfunc,
        sort_method,
        verbose,
        show_plot,
        save_dir,
        overwrite,
    )
    fig = plotter.plot()
    plotter.show_and_save(fig, "heatmap")
    return fig, plotter
