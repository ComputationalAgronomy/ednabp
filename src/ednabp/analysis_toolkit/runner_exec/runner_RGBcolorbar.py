from itertools import combinations
import numpy as np
import os
import pandas as pd
from PIL import Image
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy.spatial.distance import squareform
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale, Normalizer
from typing import override

from ..runner_build import base_runner


class RGBColorbarRunner(base_runner.AbundanceRunner):

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
        """
        Write the richness or abundance data to a CSV file.
        Richness data: number of unique species for each target per sample_id. (e.g. How many different species got in sampleA?)
        Abundance data: sum of read for each taxa per sample_id. (e.g. How many sequences got in sampleA?)

        :param write_type: The type of data to write. Can be "richness" or "abundance".
        :param taxa_level: The name of the level to target (e.g., species, family, etc.).
        :param unit_level: The name of the taxa level to group by (e.g., genes, species, haplotype, etc.). Default is "species".
        :param save_dir: The directory to save the CSV file.
        :param normalize: Whether to normalize the abundance data.
        :param sample_id_list: A list of sample IDs to write.
        """
        if write_type not in ["richness", "abundance", "betaPCA", "betaNMDS"]:
            raise ValueError("Invalid write_type. Must be 'richness', 'abundance' or 'betaPCA'.")

        os.makedirs(save_dir, exist_ok=True)

        self._load_sample_id_list(sample_id_list)

        for sample_id in self.sample_id_used:
            self._load_abundance_dict(sample_id, taxa_level, unit_level)
            if normalize:
                self._normalize_abundance_dict()
            self._abundance_dict2df(sample_id, taxa_level)
            self._update_abundance_df() # columns: target_level, Unit(species), Counts, Sample_id

        self._add_sample_info() # columns: target_level, Unit(species), Counts, Sample_id, Site, Year, Month, Sample

        self._filter_abundance_df()

        if write_type == "richness":
            self.df=self.abundance_df.groupby([taxa_level, "Site", "Year", "Month", "Sample"])["Counts"].nunique().reset_index()
            self.df.to_csv(os.path.join(save_dir, f'{unit_level}_richness.csv'), index=False)
        elif write_type == "abundance":
            self.df = self.abundance_df.groupby([taxa_level, "Site", "Year", "Month", "Sample"])["Counts"].sum().reset_index()
            self.df.to_csv(os.path.join(save_dir, f'{unit_level}_abundance.csv'), index=False)
        elif write_type == "betaPCA":
            self.df=self.abundance_df.groupby(["Site", "Year", "Month", "Sample"])[taxa_level].apply(list).reset_index(name="community")

            community_list = self.df["community"].to_list()
            distance = [ColorbarRunner._calc_sorensen(i, j) for (i, j) in combinations(community_list, 2)]
            distance_matrix = squareform(distance)

            pca = PCA(n_components=3)
            pca_axis = np.array(pca.fit_transform(distance_matrix))

            pca_rgb = ((pca_axis - pca_axis.min()) * (1/(pca_axis.max() - pca_axis.min())))
            self.df["rgb"] = [",".join(map(str, rgb)) for rgb in pca_rgb]
            self.df.to_csv(os.path.join(save_dir, f'{unit_level}_betaPCA.csv'), index=False)

            # fig = plt.figure()
            # ax = fig.add_subplot(projection='3d')
            # ax.scatter(pca_axis[:, 0], pca_axis[:, 1], pca_axis[:, 2])
            # ax.set_xlabel('PC 1')
            # ax.set_ylabel('PC 2')
            # ax.set_zlabel('PC 3')

            # plt.scatter(pca_axis[:, 0], pca_axis[:, 1])
            # plt.xlabel("Axis 1")
            # plt.ylabel("Axis 2")
            # plt.title("PCA of Beta diversity")
            # plt.show()

            # xnew = np.linspace(1, 11, 100) 
            # spl = make_interp_spline(list(range(1, 11)), pca.explained_variance_ratio_ , k=3)  # type: BSpline
            # ratio_smooth = spl(xnew)
            
            # plt.plot(xnew, ratio_smooth)
            # plt.title(f"AUC: {sum(pca.explained_variance_ratio_):.2f}")
            # plt.xlabel("Number of components")
            # plt.ylabel("Explained variance ratio")
            # plt.xticks(range(1, 11))
            # plt.show()
        else:
            self.df=self.abundance_df.groupby(["Site", "Year", "Month", "Sample"])[taxa_level].apply(list).reset_index(name="community")

            community_list = self.df["community"].to_list()
            distance = [ColorbarRunner._calc_sorensen(i, j) for (i, j) in combinations(community_list, 2)]
            distance_matrix = squareform(distance)

            nmds = manifold.MDS(n_components=3,
                                metric=False,
                                max_iter=3000,
                                eps=1e-12,
                                dissimilarity="precomputed",
                                random_state=42,
                                n_jobs=1,
                                n_init=1,)
            nmds_axis = np.array(nmds.fit_transform(distance_matrix))

            nmds_rgb = ((nmds_axis - nmds_axis.min()) * (1/(nmds_axis.max() - nmds_axis.min())))
            self.df["rgb"] = [",".join(map(str, rgb)) for rgb in nmds_rgb]
            self.df.to_csv(os.path.join(save_dir, f'{unit_level}_betaNMDS.csv'), index=False)

            # fig = plt.figure()
            # ax = fig.add_subplot(projection='3d')
            # ax.scatter(pca_axis[:, 0], pca_axis[:, 1], pca_axis[:, 2])
            # ax.set_xlabel('Axis 1')
            # ax.set_ylabel('Axis 2')
            # ax.set_zlabel('Axis 3')

            # plt.scatter(nmds_axis[:, 0], nmds_axis[:, 1])
            # plt.xlabel("Axis 1")
            # plt.ylabel("Axis 2")
            # plt.title("NMDS of Beta diversity")
            # plt.show()

        self.analysis_type = "Write species diversity to csv"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "write_type": write_type,
                "taxa_level": taxa_level,
            }
        )

    @base_runner.log_execution("Plot heatmap", "plot_heatmap.log")
    def run_plot(self,
            csv_path: str,
            x_categories: list[str],
            save_dir: str = None,
            dereplicate: bool = False,
            colorscale: str = "tempo"
        ):
        """
        Plot a heatmap of the data.

        :param csv_path: The path to the CSV file containing the data.
        :param x_categories: The categories to use for the x-axis. Can be a single string or a list of strings.
        :param save_dir: The directory to save the .HTML file. If not provided, the output will not be saved. Default is None.
        :param dereplicate: Whether to dereplicate by "Sample". Default is False.
        """
        self.df = pd.read_csv(csv_path)

        if dereplicate:
            self.df = self.df.drop("Sample", axis=1)

        self.pivot_table = pd.pivot_table(self.df,
            values="Counts",
            index=self.df.columns[0],
            columns=["Site", "Year", "Month", "Sample"],
            aggfunc="mean"
        )
        self.pivot_table = self.pivot_table.fillna(0)

        self._sort_index()

        self._load_heatmap_x(x_categories)

        self.fig = go.Figure(
            data=go.Heatmap(
                z=np.array(self.pivot_table)[self.s_index,:],
                y=self.pivot_table.index[self.s_index],
                x=self.x,
                colorscale=colorscale
            )
        )
        self.fig.update_layout(yaxis_title=None, width=1000, height=600)
        self.fig.show()

        if save_dir:
            self._save_html("Heatmap", save_dir, os.path.basename(csv_path).split(".")[0])

        self.analysis_type = "Plot heatmap"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "csv_path": csv_path,
                "x_categories": x_categories,
                "dereplicate": dereplicate,
            }
        )

    def _calc_sorensen(x, y):
        """
        Calculate the Sørensen index between two arrays.
        """
        intersection = len(set(x) & set(y))
        union = len(set(x)) + len(set(y))
        return 2 * intersection / union


