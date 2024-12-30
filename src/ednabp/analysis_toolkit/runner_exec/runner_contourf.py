from geokrige.methods import OrdinaryKriging
from geokrige.tools import TransformerGDF
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import pandas as pd

from ..runner_build import base_runner


class ContourfRunner(base_runner.AbundanceRunner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)

    def run_write(self,
            write_type: str = "abundance",
            taxa_level: str = "species",
            unit_level: str = "species",
            save_dir: str = ".",
            normalize: bool = False,
            sample_id_list: list[str] = []
        ):
        return super().run_write(
            write_type=write_type,
            taxa_level=taxa_level,
            unit_level=unit_level,
            save_dir=save_dir,
            normalize=normalize,
            sample_id_list=sample_id_list
        )

    def run_plot(self,
                 csv_path: str,
                 shp_path: str,
                 save_dir: str = None,
                 grid_density: float = 2.0,
                 value_step: float = 0.2,
                 cmap: str = "viridis"
                 ):
        self._load_csv(csv_path)

        self._model_interpolation()

        self._load_shp(shp_path)

        self._transform_grid_and_mask(grid_density)

        self._predict()

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            save_dir = os.path.join(save_dir, f"{os.path.basename(csv_path).split(".")[0]}.png")

        self._plot(value_step, cmap, save_dir)

    def _load_csv(self, csv_path: str):
        self.logger.info(f"Loading data from: {csv_path}")
        df = pd.read_csv(csv_path)
        df = df.filter(items=["Longitude", "Latitude", "Counts"])
        data = df.to_numpy()

        self.lon_lat = np.column_stack([data[:, 0], data[:, 1]]) #long/lati
        self.counts = np.array(data[:, 2])

    def _model_interpolation(self):
        self.logger.info("Fitting Kriging model...")
        self.kgn = OrdinaryKriging()
        self.kgn.load(self.lon_lat, self.counts)

        self.kgn.variogram(plot=False)
        self.kgn.fit(model='exp', plot=False)

    def _load_shp(self, shp_path: str):
        self.logger.info(f"Loading geographical map from: {shp_path}")
        self.prediction_gdf = gpd.read_file(shp_path).to_crs(crs='EPSG:4326')

    def _transform_grid_and_mask(self, grid_density: float):
        self.logger.info("Transforming geographical map to meshgrid...")
        transformer = TransformerGDF()
        transformer.load(self.prediction_gdf)

        self.grid = transformer.meshgrid(density=grid_density)
        self.mask = transformer.mask()

    def _predict(self):
        self.logger.info("Interpolating for unsampled fields...")
        self.X, self.Y = self.grid
        self.Z = self.kgn.predict(self.grid)
        self.Z[self.mask] = None

    def _plot(self, value_step: float, cmap: str, save_path: str):
        self.logger.info("Plotting...")
        fig, ax = plt.subplots()

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

        # # plot contour
        # clabels = ax.contour(self.X, self.Y, self.Z,
        #                      levels=np.arange(250, 1501, 200),
        #                      colors='k',
        #                      linewidths=0.5)
        # ax.clabel(clabels, fontsize=8)

        # add colorbar
        cax = fig.add_axes([0.93, 0.134, 0.02, 0.72])
        plt.colorbar(cbar, cax=cax, orientation='vertical')

        ax.grid(lw=0.2)
        ax.set_xlim(min(self.X[0]), max(self.X[0]))
        ax.set_ylim(min(self.Y[0]), max(self.Y[-1]))
        plt.show()
        if save_path:
            plt.savefig(save_path)
            self.logger.info(f"Saved PNG to {save_path}")