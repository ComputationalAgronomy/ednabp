import os
import warnings
from typing import TYPE_CHECKING, Annotated, Literal

import numpy as np
import plotly.graph_objects as go
import umap

from ..runner_build import DMPlotter, base_logger

if TYPE_CHECKING:
    import pandas as pd


class HeatmapPlotter(DMPlotter):
    @base_logger.prog_log("Plot heatmap")
    def plot_heatmap(
        self,
        csv_path: str,
        values: str,
        index: str,
        columns: str | Annotated[list[str], 2],
        aggfunc: Literal["mean", "sum"] = "mean",
        random_state: int = None,
        save_dir: str = None,
        overwrite: bool = False,
        **kwargs,
    ) -> "pd.DataFrame":
        """
        Plot a heatmap of the data.

        :param csv_path: The path to the CSV file containing the data.
        :param taxa_column: Column name to use for color values
        :param metric_column: Column name to use for y-axis values
        :param x_categories: The categories to use for the x-axis. Can be a single string or a list of strings.
        :param save_dir: If provided, the heatmap will be saved as a .HTML file. Default: None.
        :param overwrite: Whether to overwrite existing files. Default: False.
        """
        self._load_and_validate_data(csv_path, set(columns + [values, index]))
        self._process_data(values, index, columns, aggfunc, random_state)
        self._prepare_plot_data()
        self.kwargs = kwargs
        self._create_plot()
        self._display_and_save(save_dir, overwrite)
        return self.pivot_table

    @base_logger.prog_log("Create pivot table")
    def _process_data(self, values, index, columns, aggfunc, random_state):
        self._create_pivot_table(values, index, columns, aggfunc)
        self._sort_index(random_state)

    def _sort_index(self, random_state):
        warnings.filterwarnings("ignore", category=UserWarning, module="umap")
        self.s_index = (
            umap.UMAP(
                n_components=1, n_neighbors=15, random_state=random_state
            )
            .fit(np.array(self.pivot_table))
            .embedding_
        )
        self.s_index = np.argsort(self.s_index[:, 0])

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
        self.y = self.pivot_table.index[self.s_index]
        self.z = np.array(self.pivot_table)[self.s_index, :]

    @base_logger.prog_log("Create plot")
    def _create_plot(self):
        self.fig = go.Figure(
            data=go.Heatmap(
                z=self.z,
                y=self.y,
                x=self.x,
            )
        )

        self._update_fig_default_settings()
        self._add_fig_setting()

    def _update_fig_default_settings(self):
        FIG_DEFAULT_SETTINGS = {
            "title": None,
            "title_font_size": 24,
            "x_axis_title": None,
            "y_axis_title": None,
            "axes_title_font": 20,
            "show_xticks": True,
            "show_yticks": True,
            "axes_tick_font": 18,
            "showlegend": False,
            "legend_font": 15,
            "legend_x_position": 1.05,
            "legend_y_position": 1.0,
            "legend_orientation": "h",
            "legend_traceorder": "normal",
            "autosize": True,
            "width": None,
            "height": None,
            "colorscale": "Viridis",
            "colorbar_title": None,
            "colorbar_title_font_size": 20,
            "colorbar_tick_font_size": 16,
            "colorbar_len": 1.0,
            "colorbar_thickness": 20,
            "colorbar_x": 1.02,
            "colorbar_y": 0.5,
        }

        for key, value in self.kwargs.items():
            FIG_DEFAULT_SETTINGS[key] = value

        self.fig_sets = FIG_DEFAULT_SETTINGS

    def _add_fig_setting(self):
        self.fig.update_xaxes(
            title={
                "text": self.fig_sets["x_axis_title"],
                "font": {"size": self.fig_sets["axes_title_font"]},
            },
            tickfont={"size": self.fig_sets["axes_tick_font"]},
            showticklabels=self.fig_sets["show_xticks"],
        )

        self.fig.update_yaxes(
            title={
                "text": self.fig_sets["y_axis_title"],
                "font": {"size": self.fig_sets["axes_title_font"]},
            },
            tickfont={"size": self.fig_sets["axes_tick_font"]},
            showticklabels=self.fig_sets["show_yticks"],
        )

        self.fig.update_traces(
            colorscale=self.fig_sets["colorscale"],
            colorbar={
                "title": self.fig_sets["colorbar_title"],
                "titlefont": {
                    "size": self.fig_sets["colorbar_title_font_size"]
                },
                "tickfont": {"size": self.fig_sets["colorbar_tick_font_size"]},
                "len": self.fig_sets["colorbar_len"],
                "thickness": self.fig_sets["colorbar_thickness"],
                "x": self.fig_sets["colorbar_x"],
                "y": self.fig_sets["colorbar_y"],
            },
        )

        self.fig.update_layout(
            title={
                "text": self.fig_sets["title"],
                "font": {"size": self.fig_sets["title_font_size"]},
            }
            if self.fig_sets["title"]
            else None,
            autosize=self.fig_sets["autosize"],
            width=self.fig_sets["width"],
            height=self.fig_sets["height"],
            showlegend=self.fig_sets["showlegend"],
            legend={
                "x": self.fig_sets["legend_x_position"],
                "y": self.fig_sets["legend_y_position"],
                "traceorder": self.fig_sets["legend_traceorder"],
                "orientation": self.fig_sets["legend_orientation"],
                "font": {"size": self.fig_sets["legend_font"]},
            },
        )

    @base_logger.prog_log("Display and save (if 'save_dir' provided)")
    def _display_and_save(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            HeatmapPlotter._save_csv(self, save_dir, overwrite)
            HeatmapPlotter._save_plot(self, save_dir, overwrite)

    def _save_csv(self, save_dir, overwrite):
        csv_path = os.path.join(save_dir, "heatmap.csv")
        if os.path.exists(csv_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {csv_path}. Stop saving."
            )
            return
        self.pivot_table.to_csv(csv_path)
        self.logger.info(f"CSV saved to: {csv_path}")

    def _save_plot(self, save_html_dir: str, overwrite: bool):
        fig_path = os.path.join(save_html_dir, "heatmap.html")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {fig_path}. Stop saving."
            )
            return
        self.fig.write_html(fig_path)
        self.logger.info(f"Heatmap saved to: {fig_path}")
