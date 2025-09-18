from typing import Literal

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from . import base_plotter


class BarchartPlotter(base_plotter.Plotter):
    def __init__(
        self,
        df,
        index,
        columns,
        aggfunc,
        verbose,
        show_plot,
        save_dir,
        overwrite,
    ):
        self.read_df(df)
        self.index = index
        self.columns = columns
        self.aggfunc = aggfunc
        self.pivot_table = None
        super().__init__(verbose, show_plot, save_dir, overwrite)

    def get_pivot_table(self):
        self.pivot_table = pd.pivot_table(
            self.df,
            values=base_plotter.VALUE_COLUMN,
            index=self.index,
            columns=self.columns,
            aggfunc=self.aggfunc,
            fill_value=0,
        )

        row_sums = self.pivot_table.sum(axis=1)
        self.pivot_table = self.pivot_table.loc[
            row_sums.sort_values(ascending=False).index
        ]

    def get_x_y_color(self):
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
                self.x = [
                    [
                        "<br>".join(list(map(str, x_val))[::-1])
                        for x_val in zip(*self.x[:-1], strict=False)
                    ],
                    self.x[-1],
                ]

        self.y = np.array(self.pivot_table)
        self.color = self.pivot_table.index

    def plot(self):
        self.get_pivot_table()
        self.get_x_y_color()

        fig = go.Figure()
        for y_vals, c in zip(self.y, self.color, strict=False):
            fig.add_bar(x=self.x, y=y_vals, name=c)

        fig.update_layout(barmode="stack", showlegend=True)
        return fig


def barchart(
    df: str | pd.DataFrame,
    index: str,
    columns: str | list[str],
    aggfunc: Literal["mean", "sum"] = "mean",
    verbose: bool = False,
    show_plot: bool = True,
    save_dir: str | None = None,
    overwrite: bool = False,
) -> tuple[go.Figure, BarchartPlotter]:
    plotter = BarchartPlotter(
        df=df,
        index=index,
        columns=columns,
        aggfunc=aggfunc,
        verbose=verbose,
        show_plot=show_plot,
        save_dir=save_dir,
        overwrite=overwrite,
    )
    fig = plotter.plot()
    plotter.show_and_save(fig, "barchart")
    return fig, plotter
