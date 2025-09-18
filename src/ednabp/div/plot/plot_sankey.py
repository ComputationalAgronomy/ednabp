from collections import defaultdict
from typing import Literal

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from . import base_plotter


def human_format(num):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return "{:.2f}{}".format(num, ["", "K", "M", "G", "T", "P"][magnitude])


def hex_to_rgba(hex_color, alpha=1.0):
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


class SankeyPlotter(base_plotter.Plotter):
    def __init__(
        self,
        df,
        categories,
        aggfunc,
        verbose,
        show_plot,
        save_dir,
        overwrite,
    ):
        self.read_df(df)
        self.categories = categories
        self.aggfunc = aggfunc
        super().__init__(verbose, show_plot, save_dir, overwrite)

    def prepare_data(self):
        sankey_df = pd.DataFrame()
        for i in range(len(self.categories) - 1):
            temp_df = self.df[
                [
                    self.categories[i],
                    self.categories[i + 1],
                    base_plotter.VALUE_COLUMN,
                ]
            ]
            temp_df.columns = ["source", "target", base_plotter.VALUE_COLUMN]
            sankey_df = pd.concat([sankey_df, temp_df])

        sankey_df = (
            sankey_df.groupby(["source", "target"])
            .agg({base_plotter.VALUE_COLUMN: self.aggfunc})
            .reset_index()
        )

        self.labels = []
        self.label_categories = []
        for category in self.categories:
            uniq_labels = list(np.unique(self.df[category]))
            self.labels += uniq_labels
            self.label_categories += [category] * len(uniq_labels)

        self.source = sankey_df["source"].apply(lambda x: self.labels.index(x))
        self.target = sankey_df["target"].apply(lambda x: self.labels.index(x))
        self.count = sankey_df[base_plotter.VALUE_COLUMN]

    def create_node_labels(self):
        self.node_dict = {}

        source_values = defaultdict(int)
        for s, c in zip(self.source, self.count, strict=False):
            source_values[s] += c
        self.node_dict.update(source_values)

        target_values = defaultdict(int)
        for t, c in zip(self.target, self.count, strict=False):
            target_values[t] += c
        self.node_dict.update(target_values)

        label_values = [human_format(v) for v in self.node_dict.values()]
        self.node_label_list = [
            f"{n} - {v}"
            for n, v in zip(self.labels, label_values, strict=False)
        ]

    def get_colors(
        self,
        color_link_by,
        priority=None,
    ):
        colors = px.colors.qualitative.Plotly
        if len(colors) < len(self.labels):
            colors += px.colors.qualitative.D3
        if len(colors) < len(self.labels):
            colors += px.colors.qualitative.G10

        node_colors = [colors[i] for i in range(len(self.labels))]

        if color_link_by == "source":
            link_colors = [node_colors[idx] for idx in self.source]
        elif color_link_by == "target":
            link_colors = [node_colors[idx] for idx in self.target]
        elif color_link_by == "priority":
            link_colors = self.get_priority_colors(
                priority,
                self.source,
                self.target,
                self.label_categories,
                node_colors,
            )
        elif color_link_by == "no_color":
            link_colors = None
        else:
            raise ValueError(
                "color_link_by must be 'source', 'target', 'priority' or 'no_color'."
            )

        self.node_colors = [
            hex_to_rgba(color, alpha=1.0) for color in node_colors
        ]
        if link_colors:
            self.link_colors = [
                hex_to_rgba(color, alpha=0.5) for color in link_colors
            ]

    def get_priority_colors(self, priority, node_colors):
        if not set(priority).issubset(set(self.label_categories)):
            raise ValueError(
                f"Unknown categories: {list(set(priority) - set(self.label_categories))}"
            )

        priority_map = {category: idx for idx, category in enumerate(priority)}
        link_colors = []

        for src_idx, tgt_idx in zip(self.source, self.target, strict=False):
            src_category = self.label_categories[src_idx]
            tgt_category = self.label_categories[tgt_idx]
            src_priority = priority_map.get(src_category)
            tgt_priority = priority_map.get(tgt_category)

            if src_priority > tgt_priority:
                link_colors.append(node_colors[tgt_idx])
            elif src_priority < tgt_priority:
                link_colors.append(node_colors[src_idx])
            else:
                raise ValueError("Misconnect nodes from the same categories!")

        return link_colors

    def plot(self, color_link_by="priority", priority=None):
        self.prepare_data()
        self.create_node_labels()
        self.get_colors(color_link_by, priority)

        node_args = {
            "pad": 15,
            "thickness": 20,
            "line": {"color": "black", "width": 0.5},
            "label": self.node_label_list,
            "color": self.node_colors,
        }
        link_args = {
            "source": self.source,
            "target": self.target,
            "value": self.count,
        }
        if hasattr(self, "link_colors"):
            link_args["color"] = self.link_colors

        fig = go.Figure(
            data=[
                go.Sankey(
                    node=node_args,
                    link=link_args,
                )
            ]
        )

        for x_coordinate, column_name in enumerate(self.categories):
            fig.add_annotation(
                x=x_coordinate,
                y=1.2,
                xref="x",
                yref="paper",
                text=column_name,
                showarrow=False,
                font={"size": 20, "color": "black"},
                align="left",
            )

        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        fig.update_layout(
            autosize=True,
            paper_bgcolor="rgba(255,255,255,255)",
            plot_bgcolor="rgba(255,255,255,255)",
            xaxis={"showgrid": False},
            yaxis={"showgrid": False},
            font_color="black",
            font_size=16,
        )

        return fig


def sankey(
    df: str | pd.DataFrame,
    categories: list[str],
    aggfunc: Literal["mean", "sum"] = "mean",
    color_link_by: Literal[
        "source", "target", "priority", "no_color"
    ] = "no_color",
    priority: list[str] | None = None,
    verbose: bool = False,
    show_plot: bool = True,
    save_dir: str = None,
    overwrite: bool = False,
) -> tuple[go.Figure, SankeyPlotter]:
    plotter = SankeyPlotter(
        df=df,
        categories=categories,
        aggfunc=aggfunc,
        verbose=verbose,
        show_plot=show_plot,
        save_dir=save_dir,
        overwrite=overwrite,
    )
    fig = plotter.plot(color_link_by, priority)
    plotter.show_and_save(fig, "sankey")

    return fig, plotter
