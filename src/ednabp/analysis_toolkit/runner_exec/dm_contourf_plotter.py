import os

from geokrige.methods import OrdinaryKriging
from geokrige.tools import TransformerGDF
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pandas as pd

from ..runner_build import DMPlotter, base_logger


class ContourPlotter(DMPlotter):

    @base_logger.prog_log("Plot contour")
    def plot_contour(self,
            csv_path: str,
            shp_path: str,
            metric_column: str,
            grid_density: float = 2.0,
            value_step: float = 0.2,
            cmap: str = "viridis",
            save_dir: str = None,
            overwrite: bool = False,
        ):
        """
        Create a contour plot from geographic data with shapefile overlay.

        :param csv_path: Path to CSV file containing geographic data
        :param shp_path: Path to shapefile for map overlay
        :param metric_column: Name of the column containing values to plot
        :param grid_density: Density of the interpolation grid. Default: 2.0.
        :param value_step: Step size between contour levels. Default: 0.2.
        :param cmap: Matplotlib colormap name. Default: "viridis".
        :param save_dir: If provided, the contour will be saved as a .PNG file. Default: None.
        :param overwrite: Whether to overwrite existing files. Default: False.
        """
        self._load_and_validate_data(csv_path, metric_column)
        self._process_data(shp_path, grid_density)
        self._prepare_plot_data()
        self._create_plot(value_step, cmap)
        super()._display_and_save_plot(save_dir, overwrite)

    @base_logger.prog_log("Load and validate data")
    def _load_and_validate_data(self, csv_path: str, metric_column: str):
        self.df = pd.read_csv(csv_path)
        required_columns = ["Longitude", "Latitude", metric_column]

        missing_columns = set(required_columns) - set(self.df.columns)
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        if not self.df["Longitude"].between(-180, 180).all():
            raise ValueError("Longitude values must be between -180 and 180")
        if not self.df["Latitude"].between(-90, 90).all():
            raise ValueError("Latitude values must be between -90 and 90")
        if (self.df[required_columns[2]] < 0).any():
            raise ValueError(f"{required_columns[2]} values cannot be negative")

        data = self.df[required_columns].to_numpy()
        self.lon_lat = np.column_stack([data[:, 0], data[:, 1]]) #long/lati
        self.counts = np.array(data[:, 2])

    def _process_data(self, shp_path, grid_density):
        self._model_interpolation()
        self._load_shp(shp_path)
        self._transform_grid_and_mask(grid_density)

    @base_logger.prog_log("Model interpolation")
    def _model_interpolation(self):
        self.kgn = OrdinaryKriging()
        self.kgn.load(self.lon_lat, self.counts)

        self.kgn.variogram(plot=False)
        self.kgn.fit(model='exp', plot=False)

    @base_logger.prog_log("Load geographical map")
    def _load_shp(self, shp_path: str):
        self.prediction_gdf = gpd.read_file(shp_path).to_crs(crs='EPSG:4326')

    @base_logger.prog_log("Transform geographical map to meshgrid")
    def _transform_grid_and_mask(self, grid_density: float):
        transformer = TransformerGDF()
        transformer.load(self.prediction_gdf)

        self.grid = transformer.meshgrid(density=grid_density)
        self.mask = transformer.mask()

    @base_logger.prog_log("Interpolate for unsampled fields")
    def _prepare_plot_data(self):
        self.X, self.Y = self.grid
        self.Z = self.kgn.predict(self.grid)
        self.Z[self.mask] = None

    @base_logger.prog_log("Create plot")
    def _create_plot(self, value_step: float, cmap: str):
        self.fig, ax = plt.subplots()
        # plot geographical map
        self.prediction_gdf.plot(facecolor='none',
                                 edgecolor='black',
                                 linewidth=1.5,
                                 zorder=5,
                                 ax=ax)
        # plot contourf
        cbar = ax.contourf(self.X, self.Y, self.Z,
                           cmap=cmap,
                           levels=np.arange(0, max(self.counts), value_step),
                           extend='neither')
        # add colorbar
        cax = self.fig.add_axes([0.93, 0.134, 0.02, 0.72])
        self.fig.colorbar(cbar, cax=cax, orientation='vertical')

        ax.grid(lw=0.2)
        ax.set_xlim(min(self.X[0]), max(self.X[0]))
        ax.set_ylim(min(self.Y[0]), max(self.Y[-1]))
        # # plot contour
        # clabels = ax.contour(self.X, self.Y, self.Z,
        #                      levels=np.arange(250, 1501, 200),
        #                      colors='k',
        #                      linewidths=0.5)
        # ax.clabel(clabels, fontsize=8)

    def _save_plot(self, save_png_dir: str, overwrite: bool):
        fig_path = os.path.join(save_png_dir, "contour.png")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(f"WARNING: File already exists: {fig_path}. Stop saving.")
            return
        self.fig.savefig(fig_path)
        self.logger.info(f"Contour saved to: {fig_path}")