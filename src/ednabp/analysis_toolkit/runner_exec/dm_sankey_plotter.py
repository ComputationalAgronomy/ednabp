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
    ):
        """
        Create a rank correlation plot from a CSV file.

        :param csv_path: Path to CSV file containing data
        :param metric_column: Name of the column containing values to plot
        :param save_dir: If provided, the rank correlation plot will be saved as a .PNG file. Default: None.
        :param overwrite: Whether to overwrite existing files. Default: False.
        """
        SankeyPlotter._load_and_validate_data(
            self, csv_path, set(categories + [values])
        )
        SankeyPlotter._process_data(self, values, categories, aggfunc)
        SankeyPlotter._prepare_plot_data(self, values, categories)
        SankeyPlotter._create_plot(self, categories)
        SankeyPlotter._display_and_save(self, save_dir, overwrite)
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
        return "%.2f%s" % (num, ["", "K", "M", "G", "T", "P"][magnitude])

    @base_logger.prog_log("Create plot")
    def _create_plot(self, categories):
        self.fig = go.Figure(
            data=[
                go.Sankey(
                    node={
                        "pad": 15,
                        "thickness": 20,
                        "line": {"color": "black", "width": 0.5},
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
        # Adds 1st,2nd month on top,x_coordinate is 0 - 5 integers,column #name is specified by the list we passed
        for x_coordinate, column_name in enumerate(categories):
            self.fig.add_annotation(
                x=x_coordinate,  # Plotly recognizes 0-5 to be the x range.
                y=1.2,  # y value above 1 means above all nodes
                xref="x",
                yref="paper",
                text=column_name,
                showarrow=False,
                font={"size": 20, "color": "black"},
                align="left",
            )
        self.fig.update_xaxes(showticklabels=False)
        self.fig.update_yaxes(showticklabels=False)
        self.fig.update_layout(
            autosize=True,
            paper_bgcolor="rgba(255,255,255,255)",
            plot_bgcolor="rgba(255,255,255,255)",
            xaxis={"showgrid": False},
            yaxis={"showgrid": False},
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
