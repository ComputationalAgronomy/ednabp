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
        colorscale: str = "tempo",
        save_dir: str = None,
        overwrite: bool = False,
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
        self._create_plot(colorscale)
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
    def _create_plot(self, colorscale):
        self.fig = make_subplots(
            rows=1,
            cols=3,
            subplot_titles=[
                "Coefficient Matrix",
                "Significant Matrix",
                "Union(tril)-Intersection(triu) Matrix",
            ],
            shared_xaxes=False,
            shared_yaxes=True,
            horizontal_spacing=0.1,
        )
        self.fig.add_trace(
            go.Heatmap(
                z=self.coef_z,
                y=self.x,
                x=self.x,
                colorbar={"title": "Coef", "x": 0.29},
                colorscale=colorscale,
            ),
            row=1,
            col=1,
        )
        self.fig.add_trace(
            go.Heatmap(
                z=self.sign_z,
                y=self.x,
                x=self.x,
                colorbar={"title": "Sign", "x": 0.655},
                colorscale=colorscale,
            ),
            row=1,
            col=2,
        )
        self.fig.add_trace(
            go.Heatmap(
                z=self.union_inter_z,
                y=self.x,
                x=self.x,
                colorbar={"title": "Union-Inter"},
                colorscale=colorscale,
            ),
            row=1,
            col=3,
        )
        self.fig.update_layout(
            title_text="Comparison of Matrices",
            width=1500,
            height=500,
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
