from typing import Literal

import pandas as pd
import plotly.graph_objects as go

from . import base_plotter


class RankcorrPlotter(base_plotter.Plotter):
    def __init__(
        self,
        df,
        values,
        index,
        columns,
        aggfunc,
        rcorr,
        alpha,
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
        self.rcorr = rcorr
        self.alpha = alpha
        self.pivot_table = None
        self.corr_matrix = None
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

    def calculate_correlation(self):
        self.corr_matrix = self.pivot_table.corr(method=self.rcorr)

    def prepare_plot_data(self):
        column_names = self.corr_matrix.columns.names
        if len(column_names) == 1:
            self.labels = self.corr_matrix.columns
        else:
            self.labels = []
            for column in column_names:
                todrop_columns = column_names.copy()
                todrop_columns.remove(column)
                self.labels.append(
                    self.corr_matrix.columns.droplevel(todrop_columns)
                )
            if len(column_names) > 2:
                self.config.logger.warning(
                    "More than two-level categorical x-axis is not yet available in Plotly yet. This is a substitute implementation that combines the first n-1 categories into the first level."
                )
                self.labels = [
                    [
                        "<br>".join(list(map(str, x))[::-1])
                        for x in zip(*self.labels[:-1], strict=False)
                    ],
                    self.labels[-1],
                ]
        self.z = self.corr_matrix.values

    def plot(self):
        self.create_pivot_table()
        self.calculate_correlation()
        self.prepare_plot_data()

        fig = go.Figure(
            data=go.Heatmap(
                z=self.z,
                x=self.labels,
                y=self.labels,
                colorscale="RdBu",
                zmid=0,
            )
        )

        fig.update_layout(
            title="Rank Correlation Matrix",
            xaxis_title="Variables",
            yaxis_title="Variables",
        )

        return fig


def rankcorr(
    df: str | pd.DataFrame,
    values: str,
    index: str,
    columns: str | list[str],
    aggfunc: Literal["mean", "sum"] = "mean",
    rcorr: Literal["kendall", "spearman"] = "kendall",
    alpha: float = 0.05,
    verbose: bool = False,
    show_plot: bool = True,
    save_dir: str | None = None,
    overwrite: bool = False,
) -> tuple[go.Figure, RankcorrPlotter]:
    plotter = RankcorrPlotter(
        df,
        values,
        index,
        columns,
        aggfunc,
        rcorr,
        alpha,
        verbose,
        show_plot,
        save_dir,
        overwrite,
    )
    fig = plotter.plot()
    plotter.show_and_save(fig, "rankcorr")
    return fig, plotter
