import os
from typing import Annotated, Literal

import numpy as np
import plotly.graph_objs as go

from ..runner_build import DMPlotter, base_logger


class BarchartPlotter(DMPlotter):
    @base_logger.prog_log("Plot barchart")
    def plot_barchart(
        self,
        csv_path: str,
        values: str,
        index: str,
        columns: str | Annotated[list[str], 2],
        aggfunc: Literal["mean", "sum"] = "mean",
        save_dir: str | None = None,
        overwrite: bool = False,
        **kwargs,
    ):
        """
        Plot a barchart to visualize the abundance of a level across samples.

        :param csv_path: Path to the CSV file containing the data
        :param taxa_column: Column name to use for color values
        :param metric_column: Column name to use for y-axis values
        :param save_dir: If provided, the barchart will be saved as a .HTML file. Default is None.
        :param overwrite: Whether to overwrite existing files. Default: False.
        """
        BarchartPlotter._load_and_validate_data(
            self, csv_path, set(columns + [values, index])
        )
        BarchartPlotter._process_data(self, values, index, columns, aggfunc)
        BarchartPlotter._prepare_plot_data(self)
        BarchartPlotter._create_plot(self, kwargs)
        BarchartPlotter._display_and_save(self, save_dir, overwrite)
        return self.pivot_table

    @base_logger.prog_log("Create pivot table")
    def _process_data(self, values, index, columns, aggfunc):
        self._create_pivot_table(values, index, columns, aggfunc)
        self._sort_pivot_table()

    def _sort_pivot_table(self):
        # Sort columns by sum of values (descending)
        # column_sums = self.pivot_table.sum()
        # self.pivot_table = self.pivot_table[column_sums.sort_values(ascending=False).index]

        # Sort rows by sum of values (descending)
        row_sums = self.pivot_table.sum(axis=1)
        self.pivot_table = self.pivot_table.loc[
            row_sums.sort_values(ascending=False).index
        ]

    @base_logger.prog_log("Prepare plot data")
    def _prepare_plot_data(self):
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
                self.logger.warning(
                    "WARNING: More than two-level categorical x-axis is not yet available in Plotly yet. This is a substitute implementation that combines the first n-1 categories into the first level."
                )
                self.x = [
                    [
                        "<br>".join(list(map(str, x))[::-1])
                        for x in zip(*self.x[:-1], strict=False)
                    ],
                    self.x[-1],
                ]

        self.y = np.array(self.pivot_table)
        self.color = self.pivot_table.index

    @base_logger.prog_log("Create plot")
    def _create_plot(self, kwargs):
        self.fig = go.Figure()
        for y, c in zip(self.y, self.color, strict=False):
            self.fig.add_bar(x=self.x, y=y, name=c)
        self._update_fig_default_settings(kwargs)
        self._add_fig_setting()

    def _update_fig_default_settings(self, kwargs):
        FIG_DEFAULT_SETTINGS = {
            "x_axis_title": "Sample ID",
            "y_axis_title": "Percentage (%)",
            "axes_title_font": 20,
            "axes_tick_font": 18,
            "legend_font": 15,
            "legend_x_position": 1.05,
            "legend_y_position": 1.0,
        }
        for key, value in kwargs.items():
            if key in FIG_DEFAULT_SETTINGS:
                FIG_DEFAULT_SETTINGS[key] = value

        self.fig_sets = FIG_DEFAULT_SETTINGS

    def _add_fig_setting(
        self,
    ):
        # self.fig.update_xaxes(
        #     tickmode='linear',
        #     title=dict(
        #         text=self.fig_sets["x_axis_title"],
        #         font=dict(size=self.fig_sets["axes_title_font"])
        #         ),
        #     tickfont=dict(size=self.fig_sets["axes_tick_font"])
        # )
        # self.fig.update_yaxes(
        #     title=dict(
        #         text=self.fig_sets["y_axis_title"],
        #         font=dict(size=self.fig_sets["axes_title_font"])
        #     ),
        #     tickfont=dict(size=self.fig_sets["axes_tick_font"])
        # )
        self.fig.update_layout(
            autosize=True,
            barmode="stack",
            legend={
                "x": self.fig_sets["legend_x_position"],
                "y": self.fig_sets["legend_y_position"],
                "traceorder": "normal",
                "orientation": "h",
                "font": dict(size=self.fig_sets["legend_font"]),
            },
        )

    @base_logger.prog_log("Display and save (if 'save_dir' provided)")
    def _display_and_save(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            BarchartPlotter._save_plot(self, save_dir, overwrite)
            BarchartPlotter._save_csv(self, save_dir, overwrite)

    def _save_csv(self, save_dir, overwrite):
        csv_path = os.path.join(save_dir, "barchart.csv")
        if os.path.exists(csv_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {csv_path}. Stop saving."
            )
            return
        self.pivot_table.to_csv(csv_path)
        self.logger.info(f"CSV saved to: {csv_path}")

    def _save_plot(self, save_html_dir: str, overwrite: bool):
        fig_path = os.path.join(save_html_dir, "barchart.html")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {fig_path}. Stop saving."
            )
            return
        self.fig.write_html(fig_path)
        self.logger.info(f"Barchart saved to: {fig_path}")
