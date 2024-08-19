import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from analysis_toolkit.runner_build import base_runner
from analysis_toolkit.runner_exec import runner_umap


class HdbscanRunner(base_runner.Runner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)
        self.cluster_report = []

    def run_write(self):
        return super().run_write()

    @base_runner.log_execution("Run HDBSCAN clustering", "plot_hdbscan.log")
    def run_plot(self,
            index_file: str,
            n_unit_threshold: int,
            category: str,
            save_dir: str = '.',
            min_samples: int = 5,
            min_cluster_size: int = 5,
            cluster_selection_epsilon: float = 1.0,
            alpha: float = 1.0,
            cmap: str = "rainbow"
        ) -> None:
        """
        Cluster UMAP embeddings by HDBSCAN and plot the results based on the specified category (unit, target, or all).
        (HDBSCAN parameters reference: https://hdbscan.readthedocs.io/en/latest/parameter_selection.html)

        :param index_file: Path to the UMAP index file.
        :param n_unit_threshold: Filter out units that occur less than n times in the index DataFrame. The value should be set equal to UMAP n_neighbors.
        :param category: Column name to group the units by, restricted to 'unit', 'target', or 'all'.
        :param save_dir: Directory where the output cluster report files and plots will be saved. Default is current directory.
        :param min_samples: Provide a measure of how conservative want clustering to be for HDBSCAN. Default is 5.
        :param min_cluster_size: Minimum number of samples required to consider as a cluster for HDBSCAN. Default is 5
        :param cluster_selection_epsilon: The distance threshold that clusters below the given value are not split up any further. Default is 1.0.
        :param alpha: Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative. Default is 1.0.
        :param cmap: Colormap for the plots.
        """
        if category not in ['unit', 'target', 'all']:
            raise ValueError("Invalid category. Must be 'unit', 'target', or 'all'.")

        os.makedirs(save_dir, exist_ok=True)

        index = pd.read_csv(index_file, sep='\t')
        self.index = runner_umap.UmapRunner._filter_index_by_unit_occurrence(index, n_unit_threshold)

        self._run_hdbscan_by_category(
            category=category,
            save_dir=save_dir,
            min_samples=min_samples,
            min_cluster_size=min_cluster_size,
            cluster_selection_epsilon=cluster_selection_epsilon,
            alpha=alpha,
            cmap=cmap
        )

        self.analysis_type = "hdbscan_run"
        self.results_dir = save_dir
        self.parameters.update(
            {
            "index_file": index_file,
            "n_unit_threshold": n_unit_threshold,
            "category": category,
            "min_samples": min_samples,
            "min_cluster_size": min_cluster_size,
            "cluster_selection_epsilon": cluster_selection_epsilon,
            "alpha": alpha,
            "cmap": cmap,
            }
        )

    @staticmethod
    def _fit_hdbscan(
            points: np.ndarray,
            min_samples: int,
            min_cluster_size: int,
            cluster_selection_epsilon: float,
            alpha: float
        ) -> tuple[np.ndarray, np.ndarray, int, float]:
        """
        Fit HDBSCAN clustering to the given points and return the cluster labels, whether each point is clustered, the number of clusters, and the percentage of clustered points.

        :param points: Input data points (from UMAP embeddings) for HDBSCAN clustering.
        :param min_samples: Provide a measure of how conservative want clustering to be for HDBSCAN.
        :param min_cluster_size: Minimum number of samples required to consider as a cluster for HDBSCAN.
        :param cluster_selection_epsilon: The distance threshold that clusters below the given value are not split up any further.
        :param alpha: Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative.
        """
        labels = hdbscan.HDBSCAN(
            min_samples=min_samples,
            min_cluster_size=min_cluster_size,
            cluster_selection_epsilon=cluster_selection_epsilon,
            alpha=alpha,
            allow_single_cluster=True
        ).fit_predict(points)
        clustered = (labels >= 0)
        n_c = max(labels) + 1
        noise = sum(1 for i in labels if i < 0)
        c_p = (1 - noise / len(labels)) * 100
        return labels , clustered, n_c, round(c_p, 2) # TODO(SW): Why return 4 elements while 3 of them are never used?

    def _plot_hdbscan(self,
            points: np.ndarray,
            png_path: str,
            cmap: str,
            background: str = "white",
            width: int = 800,
            height:int = 800
        ) -> None:
        """
        Plot the HDBSCAN clustering results using matplotlib.
        """
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)
        ax.set_facecolor(background)

        point_size = 300.0 / np.sqrt(points.shape[0])

        ax.scatter(
            points[~self.clustered, 0],
            points[~self.clustered, 1],
            color=(0.5, 0.5, 0.5),
            s=point_size,
            alpha=0.5
        )
        ax.scatter(
            points[self.clustered, 0],
            points[self.clustered, 1],
            c=self.labels[self.clustered],
            s=point_size,
            cmap=cmap
        )
        ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
        ax.figure.savefig(png_path)

    def _run_hdbscan(self,
            min_samples: int,
            min_cluster_size: int,
            cluster_selection_epsilon: float,
            alpha: float,
            png_path: str,
            cmap: str
        ) -> tuple[str, str, str]:
        """
        run HDBSCAN clustering and plot the results.
        """
        points = self.subindex[["umap1", "umap2"]].to_numpy()
        self.labels, self.clustered, self.numb_clus, self.clus_perc = HdbscanRunner.fit_hdbscan(points, min_samples, min_cluster_size, cluster_selection_epsilon, alpha)
        self.numb_unit = len(self.subindex["unit"].unique())
        self._plot_hdbscan(points, png_path, cmap)
        self.logger.info(f'Saved hdbscan plot to: {png_path}')

    def _update_cluster_report(self, name: str):
        """
        Update the cluster report with the current category.
        """
        self.cluster_report.append([name, self.numb_unit, self.numb_clus, self.clus_perc])

    def _write_cluster_report(self, cluster_report_path) -> None:
        """
        Write the cluster report to a TSV file.
        Columns: name, unit_counts, cluster_counts, clustered_ratio
        """
        self.cluster_report = pd.DataFrame(
            self.cluster_report,
            columns = [
                "name",
                "unit_counts",
                "cluster_counts",
                "clustered_ratio"
            ]
        )
        self.cluster_report.to_csv(cluster_report_path, sep='\t', index=False)
        self.logger.info(f'Saved cluster report to: {cluster_report_path}')

    def _run_hdbscan_by_category(self,
            category: str,
            save_dir: str,
            min_samples: int,
            min_cluster_size: int,
            cluster_selection_epsilon: float,
            alpha: float,
            cmap: str,
        ) -> None:
        """
        Run HDBSCAN clustering for a specified category and plot the results, grouped by the specified category.
        Writes a cluster report containing the number of units, clusters, and clustered ratio for each group.

        :param index: Input index DataFrame.
        :param min_samples: Provide a measure of how conservative want clustering to be for HDBSCAN.
        :param min_cluster_size: Minimum number of samples required to consider as a cluster for HDBSCAN.
        :param cluster_selection_epsilon: The distance threshold that clusters below the given value are not split up any further.
        :param alpha: Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative.
        :param category: Column name to group the units by, restricted to 'unit', 'target', or 'all'.
        :param prefix: Prefix for the output file names.
        :param save_dir: Directory to save the output files.
        :param cmap: Color map for the plot.
        """
        cluster_report_path = os.path.join(save_dir, f"{category}_cluster_report.tsv")

        if category == 'all':
            self.logger.info("Clustering for all units...")
            png_path = os.path.join(save_dir, f"all_hdbscan.png")
            self.subindex = self.index.copy()

            self._run_hdbscan(
                min_samples,
                min_cluster_size,
                cluster_selection_epsilon,
                alpha,
                png_path,
                cmap
            )

            self._update_cluster_report(category)
            self._write_cluster_report(cluster_report_path)

            return

        unique_values = np.unique(self.index[category])
        for value in unique_values:
            self.logger.info(f"Clustering for {category}: {value}...")
            png_path = os.path.join(save_dir, f"{value}_hdbscan.png")

            self.subindex = self.index[self.index[category] == value]

            self._run_hdbscan(
                min_samples,
                min_cluster_size,
                cluster_selection_epsilon,
                alpha,
                png_path,
                cmap
            )
            self._update_cluster_report(value)

        self._write_cluster_report(cluster_report_path)