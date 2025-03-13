from typing import override

import matplotlib.pyplot as plt
import numpy as np

from ..runner_build import Clusterer

# hdbscan settings (https://scikit-learn.org/1.5/modules/generated/sklearn.cluster.HDBSCAN.html)
HDBSCAN_DEFAULT_SETTINGS = {
    "min_cluster_size": 5,
    "min_samples": None,
    "metric": "euclidean",
    "p": None,
    "alpha": 1.0,
    "cluster_selection_epsilon": 0.0,
    "algorithm": "best",
    "leaf_size": 40,
    "approx_min_span_tree": True,
    "gen_min_span_tree": False,
    "core_dist_n_jobs": 4,
    "cluster_selection_method": "eom",
    "allow_single_cluster": False,
    "prediction_data": False,
    "match_reference_implementation": False,
}


class HdbClusterer(Clusterer):
    def __init__(self):
        super().__init__()

    @override
    def run(
        self,
        index_path: str | None = None,
        points: np.ndarray | None = None,
        true_labels: np.ndarray | None = None,
        **settings,
    ) -> tuple[int, int, float, float, float]:
        """
        :param index_path: Path to the index file containing the UMAP embeddings to cluster.
        :param settings: Keyword arguments for both plotting and HDBSCAN configuration:
            Metrics log settings:
                - print_metrics_log: Whether to print metrics. Default to True.
                - metrics_log_path: Path to save the calculated metrics log. Default to None (don't save).
            Plotting settings:
                - cmap (str): Colormap for plotting. Default to "Spectral".
                - background (str): Background color for plot. Default to "white".
                - width (int): Width of the output plot in pixels. Default to 800.
                - height (int): Height of the output plot in pixels. Default to 800.
                - show_plot (bool): Whether to display the plot. Default to True.
                - plot_path (str): Path to save the plot. Default to None (don't save).
            HDBSCAN settings:
                - Any parameter accepted by HDBSCAN (https://scikit-learn.org/1.5/modules/generated/sklearn.cluster.HDBSCAN.html)
        :return: A tuple containing five clustering metrics:
            - actual_num (int): The actual number of clusters set
            - cluster_num (int): The number of clusters specified
            - cluster_perc (float): The percentage of data points assigned to clusters
            - silhouette_avg (float): The average silhouette score for the clustering
            - ari (float): The Adjusted Rand Index score
        """
        self._add_default_settings(settings=settings)
        return super().run(
            index_path=index_path, points=points, true_labels=true_labels
        )

    @override
    def _add_default_settings(self, settings):
        DEFAULT_SETTINGS = HDBSCAN_DEFAULT_SETTINGS
        for key, value in DEFAULT_SETTINGS.items():
            if key not in settings:
                settings[key] = value
        self.settings = settings
        super()._add_default_settings()

    @override
    def _run_clustering(self) -> None:
        self.init_clusterer()
        self.fit_hdbscan()
        if (
            self.settings["plot_path"] is not None
            or self.settings["show_plot"] is True
        ):
            self.output_hdbscan()

    def init_clusterer(self):
        import hdbscan

        self.clusterer = hdbscan.HDBSCAN(
            min_cluster_size=self.settings["min_cluster_size"],
            min_samples=self.settings["min_samples"],
            metric=self.settings["metric"],
            p=self.settings["p"],
            alpha=self.settings["alpha"],
            cluster_selection_epsilon=self.settings[
                "cluster_selection_epsilon"
            ],
            algorithm=self.settings["algorithm"],
            leaf_size=self.settings["leaf_size"],
            approx_min_span_tree=self.settings["approx_min_span_tree"],
            gen_min_span_tree=self.settings["gen_min_span_tree"],
            core_dist_n_jobs=self.settings["core_dist_n_jobs"],
            cluster_selection_method=self.settings["cluster_selection_method"],
            allow_single_cluster=self.settings["allow_single_cluster"],
            prediction_data=self.settings["prediction_data"],
            match_reference_implementation=self.settings[
                "match_reference_implementation"
            ],
        )

    def fit_hdbscan(self) -> None:
        self.cluster_labels = self.clusterer.fit_predict(self.points)
        self.clustered = self.cluster_labels >= 0

    def output_hdbscan(self) -> None:
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(
            figsize=(
                self.settings["width"] / dpi,
                self.settings["height"] / dpi,
            )
        )
        ax = fig.add_subplot(111)
        ax.set_facecolor(self.settings["background"])

        point_size = 300.0 / np.sqrt(self.points.shape[0])

        ax.scatter(
            self.points[~self.clustered, 0],
            self.points[~self.clustered, 1],
            color=(0.5, 0.5, 0.5),
            s=point_size,
            alpha=0.5,
        )
        ax.scatter(
            self.points[self.clustered, 0],
            self.points[self.clustered, 1],
            c=self.cluster_labels[self.clustered],
            s=point_size,
            cmap=self.settings["cmap"],
        )
        ax.tick_params(
            bottom=False, left=False, labelbottom=False, labelleft=False
        )
        if self.settings["plot_path"] is not None:
            ax.figure.savefig(self.settings["plot_path"])
        if self.settings["show_plot"] is True:
            plt.show()
