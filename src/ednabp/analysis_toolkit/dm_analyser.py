from typing import TYPE_CHECKING, Annotated, Literal, Union

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

    from .mito_data import MitoData
    from .sample_data import SampleData


class DMAnalyser:
    """
    A class integrates several data writing and plotting functionalities for performing diversity metrics analysis and visualization.
    It is designed to handle and analyze diversity metrics from the class SampleData or MitoData.

    :attributes sample_data: SampleData or MitoData (optional, only be used for data writing). The input data used for analysis and visualization.
    :attributes verbose: Whether to print verbose output during processing.
    :attributes dm_writer: A data writing class for outputting diversity metrics tables.
    :attributes barchart_plotter: A plotter class for generating bar charts.
    :attributes heatmap_plotter: A plotter class for generating heatmaps.
    :attributes rankcorr_plotter: A plotter class for generating rank correlation plots.
    :attributes sankey_plotter: A plotter class for generating Sankey diagrams.
    :attributes contour_plotter: A plotter class for generating contour plots.
    """

    def __init__(
        self,
        sample_data: Union["SampleData", "MitoData", None] = None,
        verbose: bool = True,
    ):
        self.sample_data = sample_data
        self.verbose = verbose

    def _import_dmwriter(self):
        from .runner_exec.dm_table_writer import DMWriter

        self.dm_writer = DMWriter(self.sample_data, self.verbose)

    def write_richness_table(
        self,
        save_dir: str,
        taxa_level: str,
        unit_level: str = "species",
        sample_id_list: list[str] | None = None,
        overwrite: bool = False,
    ) -> "pd.DataFrame":
        if not hasattr(self, "dm_writer"):
            self._import_dmwriter()
        return self.dm_writer.write_richness_table(
            save_dir, taxa_level, unit_level, sample_id_list, overwrite
        )

    def write_abundance_table(
        self,
        save_dir: str,
        taxa_level: str,
        process: Literal["norm", "log"] | None = None,
        sample_id_list: list[str] | None = None,
        overwrite: bool = False,
    ) -> "pd.DataFrame":
        if not hasattr(self, "dm_writer"):
            self._import_dmwriter()
        return self.dm_writer.write_abundance_table(
            save_dir, taxa_level, process, sample_id_list, overwrite
        )

    def write_detectprob_table(
        self,
        save_dir: str,
        taxa_level: str,
        detectprob_column: str | list[str] = "Sample",
        sample_id_list: list[str] | None = None,
        overwrite: bool = False,
    ) -> "pd.DataFrame":
        if not hasattr(self, "dm_writer"):
            self._import_dmwriter()
        return self.dm_writer.write_detectprob_table(
            save_dir, taxa_level, detectprob_column, sample_id_list, overwrite
        )

    def _import_barchartplotter(self):
        from .runner_exec.dm_barchart_plotter import BarchartPlotter

        self.barchart_plotter = BarchartPlotter(self.verbose)

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
    ) -> "pd.DataFrame":
        if not hasattr(self, "barchart_plotter"):
            self._import_barchartplotter()
        return self.barchart_plotter.plot_barchart(
            csv_path,
            values,
            index,
            columns,
            aggfunc,
            save_dir,
            overwrite,
            **kwargs,
        )

    def _import_heatmapplotter(self):
        from .runner_exec.dm_heatmap_plotter import HeatmapPlotter

        self.heatmap_plotter = HeatmapPlotter(self.verbose)

    def plot_heatmap(
        self,
        csv_path: str,
        values: str,
        index: str,
        columns: str | Annotated[list[str], 2],
        aggfunc: Literal["mean", "sum"] = "mean",
        colorscale: str = "tempo",
        save_dir: str = None,
        overwrite: bool = False,
    ) -> "pd.DataFrame":
        if not hasattr(self, "heatmap_plotter"):
            self._import_heatmapplotter()

        return self.heatmap_plotter.plot_heatmap(
            csv_path,
            values,
            index,
            columns,
            aggfunc,
            colorscale,
            save_dir,
            overwrite,
        )

    def _import_rankcorrplotter(self):
        from .runner_exec.dm_rankcorr_plotter import RankCorrPlotter

        self.rankcorr_plotter = RankCorrPlotter(self.verbose)

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
    ) -> "pd.DataFrame":
        if not hasattr(self, "rankcorr_plotter"):
            self._import_rankcorrplotter()
        return self.rankcorr_plotter.plot_rankcorr(
            csv_path,
            values,
            index,
            columns,
            aggfunc,
            rcorr,
            alpha,
            colorscale,
            save_dir,
            overwrite,
        )

    def _import_sankeyplotter(self):
        from .runner_exec.dm_sankey_plotter import SankeyPlotter

        self.sankey_plotter = SankeyPlotter(self.verbose)

    def plot_sankey(
        self,
        csv_path: str,
        values: str,
        categories: list[str],
        aggfunc: Literal["mean", "sum"] = "sum",
        save_dir: str = None,
        overwrite: bool = False,
    ) -> "pd.DataFrame":
        if not hasattr(self, "sankey_plotter"):
            self._import_sankeyplotter()
        return self.sankey_plotter.plot_sankey(
            csv_path, values, categories, aggfunc, save_dir, overwrite
        )

    def _import_contourplotter(self):
        from .runner_exec.dm_contourf_plotter import ContourPlotter

        self.contour_plotter = ContourPlotter(self.verbose)

    def plot_contour(
        self,
        csv_path: str,
        shp_path: str,
        metric_column: str,
        grid_density: float = 2.0,
        value_step: float = 0.2,
        cmap: str = "viridis",
        save_dir: str = None,
        overwrite: bool = False,
    ) -> "np.array":
        if not hasattr(self, "contour_plotter"):
            self._import_contourplotter()
        return self.contour_plotter.plot_contour(
            csv_path,
            shp_path,
            metric_column,
            grid_density,
            value_step,
            cmap,
            save_dir,
            overwrite,
        )
