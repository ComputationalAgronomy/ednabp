import os
from collections import defaultdict
from itertools import product
from typing import Annotated, Literal

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..runner_build import DMPlotter, base_logger


class SankeyPlotter(DMPlotter):
    @base_logger.prog_log("Create sankey diagram")
    def plot_sankey(
        self,
        csv_path: str,
        values: str,
        categories: list[str],
        aggfunc: Literal["mean", "sum"] = "sum",
        save_dir: str = None,
        overwrite: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Create a rank correlation plot from a CSV file.

        :param csv_path: Path to CSV file containing data
        :param metric_column: Name of the column containing values to plot
        :param save_dir: If provided, the rank correlation plot will be saved as a .PNG file. Default: None.
        :param overwrite: Whether to overwrite existing files. Default: False.
        """
        self._load_and_validate_data(csv_path, set(categories + [values]))
        self._process_data(values, categories, aggfunc)
        self._prepare_plot_data(values, categories)
        self.kwargs = kwargs
        self._create_plot(categories)
        self._display_and_save(save_dir, overwrite)
        return self.sankey_df

    def _process_data(self, values, categories, aggfunc):
        self.sankey_df = pd.DataFrame()
        for i in range(len(categories) - 1):
            tempDf = self.df[[categories[i], categories[i + 1], values]]
            tempDf.columns = ["source", "target", values]
            self.sankey_df = pd.concat([self.sankey_df, tempDf])
        self.sankey_df = (
            self.sankey_df.groupby(["source", "target"])
            .agg({values: aggfunc})
            .reset_index()
        )

    @base_logger.prog_log("Prepare plot data")
    def _prepare_plot_data(self, values, categories):
        node_names = list(np.unique(self.df[categories].values))
        self.source = self.sankey_df["source"].apply(
            lambda x: node_names.index(x)
        )
        self.target = self.sankey_df["target"].apply(
            lambda x: node_names.index(x)
        )
        self.count = self.sankey_df[values]

        node_dict = {}
        source_values = defaultdict(int)
        for s, c in zip(self.source, self.count, strict=False):
            source_values[s] += c
        node_dict.update(source_values)
        target_values = defaultdict(int)
        for t, c in zip(self.target, self.count, strict=False):
            target_values[t] += c
        node_dict.update(target_values)

        node_values = [
            self._human_format(node_dict[i]) for i in range(len(node_dict))
        ]

        self.label_list = [
            f"{n} - {v}" for n, v in zip(node_names, node_values, strict=False)
        ]

    @staticmethod
    def _human_format(num):
        magnitude = 0
        while abs(num) >= 1000:
            magnitude += 1
            num /= 1000.0
        # add more suffixes if you need them
        return "{:.2f}{}".format(num, ["", "K", "M", "G", "T", "P"][magnitude])

    @base_logger.prog_log("Create plot")
    def _create_plot(self, categories):
        self.fig = go.Figure(
            data=[
                go.Sankey(
                    node={
                        "label": self.label_list,
                    },
                    link={
                        "source": self.source,
                        "target": self.target,
                        "value": self.count,
                    },
                )
            ]
        )

        self._update_fig_default_settings()
        self._add_fig_setting(categories)

    def _update_fig_default_settings(self):
        FIG_DEFAULT_SETTINGS = {
            "title": None,
            "title_font_size": 24,
            "x_axis_title": None,
            "y_axis_title": None,
            "axes_title_font": 20,
            "show_xticks": False,
            "show_yticks": False,
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
            # Sankey specific settings
            "node_pad": 15,
            "node_thickness": 20,
            "node_line_color": "black",
            "node_line_width": 0.5,
            "paper_bgcolor": "rgba(255,255,255,255)",
            "plot_bgcolor": "rgba(255,255,255,255)",
            "annotation_font_size": 20,
            "annotation_font_color": "black",
            "annotation_y_position": 1.2,
        }
        for key, value in self.kwargs.items():
            FIG_DEFAULT_SETTINGS[key] = value

        self.fig_sets = FIG_DEFAULT_SETTINGS

    def _add_fig_setting(self, categories):
        self.fig.update_traces(
            node={
                "pad": self.fig_sets["node_pad"],
                "thickness": self.fig_sets["node_thickness"],
                "line": {
                    "color": self.fig_sets["node_line_color"],
                    "width": self.fig_sets["node_line_width"],
                },
            },
            selector={"type": "sankey"},
        )

        for x_coordinate, column_name in enumerate(categories):
            self.fig.add_annotation(
                x=x_coordinate,  # Plotly recognizes 0-5 to be the x range.
                y=self.fig_sets[
                    "annotation_y_position"
                ],  # y value above 1 means above all nodes
                xref="x",
                yref="paper",
                text=column_name,
                showarrow=False,
                font={
                    "size": self.fig_sets["annotation_font_size"],
                    "color": self.fig_sets["annotation_font_color"],
                },
                align="left",
            )

        self.fig.update_xaxes(
            showticklabels=self.fig_sets["show_xticks"],
            title={
                "text": self.fig_sets["x_axis_title"],
                "font": {"size": self.fig_sets["axes_title_font"]},
            },
            tickfont={"size": self.fig_sets["axes_tick_font"]},
            showgrid=False,
        )
        self.fig.update_yaxes(
            showticklabels=self.fig_sets["show_yticks"],
            title={
                "text": self.fig_sets["y_axis_title"],
                "font": {"size": self.fig_sets["axes_title_font"]},
            },
            tickfont={"size": self.fig_sets["axes_tick_font"]},
            showgrid=False,
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
            paper_bgcolor=self.fig_sets["paper_bgcolor"],
            plot_bgcolor=self.fig_sets["plot_bgcolor"],
            font_color="black",
            font_size=16,
        )

    @base_logger.prog_log("Display and save (if 'save_dir' provided)")
    def _display_and_save(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            SankeyPlotter._save_csv(self, save_dir, overwrite)
            SankeyPlotter._save_plot(self, save_dir, overwrite)

    def _save_csv(self, save_dir, overwrite):
        csv_path = os.path.join(save_dir, "sankey.csv")
        if os.path.exists(csv_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {csv_path}. Stop saving."
            )
            return
        self.sankey_df.to_csv(csv_path)
        self.logger.info(f"CSV saved to: {csv_path}")

    def _save_plot(self, save_html_dir: str, overwrite: bool):
        fig_path = os.path.join(save_html_dir, "sankey.html")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {fig_path}. Stop saving."
            )
            return
        self.fig.write_html(fig_path)
        self.logger.info(f"Heatmap saved to: {fig_path}")
