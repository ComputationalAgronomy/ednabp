import os
from itertools import product
from typing import Annotated, Literal

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import kendalltau, spearmanr

from ..runner_build import DMPlotter, base_logger


class RankCorrPlotter(DMPlotter):
    @base_logger.prog_log("Create rank correlation plot")
    def plot_rankcorr(
        self,
        csv_path: str,
        values: str,
        index: str,
        columns: str | Annotated[list[str], 2],
        aggfunc: Literal["mean", "sum"] = "mean",
        rcorr: Literal["kendall", "spearman"] = "kendall",
        alpha: float = 0.05,
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
        self._load_and_validate_data(csv_path, set(columns + [values, index]))
        self._process_data(values, index, columns, aggfunc, rcorr, alpha)
        self._prepare_plot_data()
        self.kwargs = kwargs
        self._create_plot()
        self._display_and_save(save_dir, overwrite)
        return self.corr_df

    @staticmethod
    def _calc_kendalltau(data1, data2):
        coef, p = kendalltau(data1, data2)
        return coef, p

    @staticmethod
    def _calc_spearmanr(data1, data2):
        coef, p = spearmanr(data1, data2)
        return coef, p

    @staticmethod
    def _calc_union_inter_num(data1, data2):
        data1_exist_set = {i for i, e in enumerate(data1) if e != 0}
        data2_exist_set = {i for i, e in enumerate(data2) if e != 0}
        union_num = len(data1_exist_set.union(data2_exist_set))
        inter_num = len(data1_exist_set.intersection(data2_exist_set))
        return union_num, inter_num

    @base_logger.prog_log("Create Rank Correlation Data")
    def _process_data(self, values, index, columns, aggfunc, rcorr, alpha):
        self._create_pivot_table(values, index, columns, aggfunc)
        self._process_corr_z(rcorr, alpha)
        self._create_corr_df(rcorr, alpha)

    def _process_corr_z(self, rcorr, alpha):
        func = (
            self._calc_kendalltau
            if rcorr == "kendall"
            else self._calc_spearmanr
        )
        metric_list = [
            self.pivot_table[col].to_list() for col in self.pivot_table
        ]
        self.rank_corr = [
            func(data1, data2)
            for (data1, data2) in product(metric_list, metric_list)
        ]
        self.coef_list = [v[0] for v in self.rank_corr]
        self.p_list = [v[1] for v in self.rank_corr]
        self.sign_list = [1 if v[1] < alpha else 0 for v in self.rank_corr]
        self.union_inter_num = [
            self._calc_union_inter_num(data1, data2)
            for (data1, data2) in product(metric_list, metric_list)
        ]
        self.union_inter_list = [
            self.union_inter_num[i][0] if j < k else self.union_inter_num[i][1]
            for i, (j, k) in enumerate(
                product(range(len(metric_list)), range(len(metric_list)))
            )
        ]

    def _create_corr_df(self, rcorr, alpha):
        x = [
            x
            for (x, _) in product(
                self.pivot_table.columns, self.pivot_table.columns
            )
        ]
        y = [
            y
            for (_, y) in product(
                self.pivot_table.columns, self.pivot_table.columns
            )
        ]
        self.corr_df = pd.DataFrame(
            {
                "x": x,
                "y": y,
                rcorr: self.coef_list,
                "p-value": self.p_list,
                f"significance({alpha})": self.sign_list,
                "(union, intersection)": self.union_inter_num,
            }
        )

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

        shp = int(len(self.rank_corr) ** (1 / 2))
        self.coef_z = np.array(self.coef_list).reshape(shp, shp)
        self.sign_z = np.array(self.sign_list).reshape(shp, shp)
        self.union_inter_z = np.array(self.union_inter_list).reshape(shp, shp)

    @base_logger.prog_log("Create plot")
    def _create_plot(self):
        self.fig = make_subplots(
            rows=1,
            cols=3,
            horizontal_spacing=0.1,
        )
        self.fig.add_trace(
            go.Heatmap(z=self.coef_z, y=self.x, x=self.x),
            row=1,
            col=1,
        )
        self.fig.add_trace(
            go.Heatmap(z=self.sign_z, y=self.x, x=self.x),
            row=1,
            col=2,
        )
        self.fig.add_trace(
            go.Heatmap(z=self.union_inter_z, y=self.x, x=self.x),
            row=1,
            col=3,
        )

        self._update_fig_default_settings()
        self._add_fig_setting()

    def _update_fig_default_settings(self):
        FIG_DEFAULT_SETTINGS = {
            "title": "Comparison of Matrices",
            "title_font_size": 24,
            "title_x": 0.5,
            "title_y": 0.95,
            "title_xanchor": "center",
            "title_yanchor": "top",
            "x_axis_title": None,
            "y_axis_title": None,
            "axes_title_font": 20,
            "show_xticks": True,
            "show_yticks": True,
            "show_yticks_shared": False,
            "axes_tick_font": 18,
            "showlegend": False,
            "legend_font": 15,
            "legend_x_position": 1.05,
            "legend_y_position": 1.0,
            "legend_orientation": "h",
            "legend_traceorder": "normal",
            "autosize": True,
            "width": 1500,
            "height": 500,
            # Margin settings
            "margin_left": 80,
            "margin_right": 80,
            "margin_top": 100,
            "margin_bottom": 80,
            # Background settings
            "paper_bgcolor": "white",
            "plot_bgcolor": "white",
            "font_family": "Arial",
            # Grid settings
            "xgap": 0.1,
            "ygap": 0.1,
            # Subplot title font size
            "subplot_title_font_size": 18,
            # Subplot specific settings
            "shared_xaxes": False,
            "shared_yaxes": True,
            "subplot_title_1": "Coefficient Matrix",
            "subplot_title_2": "Significant Matrix",
            "subplot_title_3": "Union(tril)-Intersection(triu) Matrix",
            # Tick modes for axes
            "tickmode_x_1": "auto",
            "tickmode_x_2": "auto",
            "tickmode_x_3": "auto",
            "tickmode_y_1": "auto",
            "tickmode_y_2": "auto",
            "tickmode_y_3": "auto",
            # Axis titles for subplots
            "x_axis_title_1": None,
            "x_axis_title_2": None,
            "x_axis_title_3": None,
            "y_axis_title_1": None,
            "y_axis_title_2": None,
            "y_axis_title_3": None,
            # Colorbar settings
            "colorbar_title_1": "Coef",
            "colorbar_x_1": 0.29,
            "colorbar_title_2": "Sign",
            "colorbar_x_2": 0.655,
            "colorbar_title_3": "Union-Inter",
            "colorbar_x_3": 0.99,
            "colorbar_thickness_1": 20,
            "colorbar_thickness_2": 20,
            "colorbar_thickness_3": 20,
            "colorbar_len_1": 0.8,
            "colorbar_len_2": 0.8,
            "colorbar_len_3": 0.8,
            # Z range settings
            "zmin_1": None,
            "zmax_1": None,
            "zmin_2": None,
            "zmax_2": None,
            "zmin_3": None,
            "zmax_3": None,
            # Colorscale options
            "colorscale_1": "RdBu",
            "colorscale_2": "RdBu",
            "colorscale_3": "RdBu",
        }

        # Update with user-provided settings
        if hasattr(self, "kwargs"):
            for key, value in self.kwargs.items():
                FIG_DEFAULT_SETTINGS[key] = value

        self.fig_sets = FIG_DEFAULT_SETTINGS

    def _add_fig_setting(self):
        self.fig.update_layout(
            grid={
                "rows": 1,
                "columns": 3,
                "pattern": "independent",
                "xgap": self.fig_sets["xgap"],
                "ygap": self.fig_sets["ygap"],
            }
        )

        self.fig.update_layout(
            annotations=[
                {
                    "text": self.fig_sets["subplot_title_1"],
                    "font": {"size": self.fig_sets["subplot_title_font_size"]},
                    "xref": "x domain",
                    "yref": "y domain",
                    "x": 0.5,
                    "y": 1.05,
                    "xanchor": "center",
                    "yanchor": "bottom",
                    "showarrow": False,
                },
                {
                    "text": self.fig_sets["subplot_title_2"],
                    "font": {"size": self.fig_sets["subplot_title_font_size"]},
                    "xref": "x2 domain",
                    "yref": "y2 domain",
                    "x": 0.5,
                    "y": 1.05,
                    "xanchor": "center",
                    "yanchor": "bottom",
                    "showarrow": False,
                },
                {
                    "text": self.fig_sets["subplot_title_3"],
                    "font": {"size": self.fig_sets["subplot_title_font_size"]},
                    "xref": "x3 domain",
                    "yref": "y3 domain",
                    "x": 0.5,
                    "y": 1.05,
                    "xanchor": "center",
                    "yanchor": "bottom",
                    "showarrow": False,
                },
            ]
        )

        # Configure shared axes
        if self.fig_sets["shared_xaxes"]:
            self.fig.update_layout(
                xaxis2={"matches": "x"}, xaxis3={"matches": "x"}
            )
        if self.fig_sets["shared_yaxes"]:
            self.fig.update_layout(
                yaxis2={"matches": "y"}, yaxis3={"matches": "y"}
            )

        for i, trace_idx in enumerate([0, 1, 2], 1):
            colorbar_settings = {
                "title": self.fig_sets[f"colorbar_title_{i}"],
                "thickness": self.fig_sets[f"colorbar_thickness_{i}"],
                "len": self.fig_sets[f"colorbar_len_{i}"],
                "x": self.fig_sets[f"colorbar_x_{i}"],
            }

            self.fig.data[trace_idx].update(
                colorbar=colorbar_settings,
                colorscale=self.fig_sets[f"colorscale_{i}"],
                zmin=self.fig_sets[f"zmin_{i}"],
                zmax=self.fig_sets[f"zmax_{i}"],
            )

        # Update x-axes for all subplots
        for i in range(1, 4):
            x_axis_name = "xaxis" if i == 1 else f"xaxis{i}"
            self.fig.update_layout(
                **{
                    x_axis_name: {
                        "title": {
                            "text": self.fig_sets[f"x_axis_title_{i}"]
                            if self.fig_sets[f"x_axis_title_{i}"] is not None
                            else self.fig_sets["x_axis_title"],
                            "font": {"size": self.fig_sets["axes_title_font"]},
                        },
                        "tickfont": {"size": self.fig_sets["axes_tick_font"]},
                        "showticklabels": self.fig_sets["show_xticks"],
                        "tickmode": self.fig_sets[f"tickmode_x_{i}"],
                    }
                }
            )

        # Update y-axes for all subplots
        for i in range(1, 4):
            y_axis_name = "yaxis" if i == 1 else f"yaxis{i}"
            show_yticks = (
                self.fig_sets["show_yticks"]
                if i == 1
                else self.fig_sets["show_yticks_shared"]
            )

            self.fig.update_layout(
                **{
                    y_axis_name: {
                        "title": {
                            "text": self.fig_sets[f"y_axis_title_{i}"]
                            if self.fig_sets[f"y_axis_title_{i}"] is not None
                            else self.fig_sets["y_axis_title"],
                            "font": {"size": self.fig_sets["axes_title_font"]},
                        },
                        "tickfont": {"size": self.fig_sets["axes_tick_font"]},
                        "showticklabels": show_yticks,
                        "tickmode": self.fig_sets[f"tickmode_y_{i}"],
                    }
                }
            )

        self.fig.update_layout(
            title={
                "text": self.fig_sets["title"],
                "font": {"size": self.fig_sets["title_font_size"]},
                "x": self.fig_sets["title_x"],
                "y": self.fig_sets["title_y"],
                "xanchor": self.fig_sets["title_xanchor"],
                "yanchor": self.fig_sets["title_yanchor"],
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
            margin={
                "l": self.fig_sets["margin_left"],
                "r": self.fig_sets["margin_right"],
                "t": self.fig_sets["margin_top"],
                "b": self.fig_sets["margin_bottom"],
            },
            paper_bgcolor=self.fig_sets["paper_bgcolor"],
            plot_bgcolor=self.fig_sets["plot_bgcolor"],
            font_family=self.fig_sets["font_family"],
        )

    @base_logger.prog_log("Display and save (if 'save_dir' provided)")
    def _display_and_save(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            RankCorrPlotter._save_csv(self, save_dir, overwrite)
            RankCorrPlotter._save_plot(self, save_dir, overwrite)

    def _save_csv(self, save_dir, overwrite):
        csv_path = os.path.join(save_dir, "rank_correlation.csv")
        if os.path.exists(csv_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {csv_path}. Stop saving."
            )
            return
        self.corr_df.to_csv(csv_path)
        self.logger.info(f"CSV saved to: {csv_path}")

    def _save_plot(self, save_html_dir: str, overwrite: bool):
        fig_path = os.path.join(save_html_dir, "rank_correlation.html")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(
                f"WARNING: File already exists: {fig_path}. Stop saving."
            )
            return
        self.fig.write_html(fig_path)
        self.logger.info(f"Heatmap saved to: {fig_path}")
