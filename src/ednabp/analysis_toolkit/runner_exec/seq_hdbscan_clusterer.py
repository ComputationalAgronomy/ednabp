from typing import override

import hdbscan
import matplotlib.pyplot as plt
import numpy as np

from ..runner_build import Clusterer

# hdbscan settings (https://scikit-learn.org/1.5/modules/generated/sklearn.cluster.HDBSCAN.html)
HDBSCAN_DEFAULT_SETTINGS = {
    "min_cluster_size": 5,
    "min_samples": None,
    "metric": 'euclidean',
    "p": None,
    "alpha": 1.0,
    "cluster_selection_epsilon": 0.0,
    "algorithm": 'best',
    "leaf_size": 40,
    "approx_min_span_tree": True,
    "gen_min_span_tree": False,
    "core_dist_n_jobs": 4,
    "cluster_selection_method": 'eom',
    "allow_single_cluster": False,
    "prediction_data": False,
    "match_reference_implementation": False,
}
class HdbClusterer(Clusterer):

    def __init__(
            self,
            index,
            **kwargs
        ):
        HdbClusterer.add_default_settings(self, settings=kwargs)
        super().__init__(index=index)

    @override
    def add_default_settings(self, settings):
        DEFAULT_SETTINGS = {
            # output settings
            "show_plot": True,
            "plot_path": None
        }
        DEFAULT_SETTINGS.update(HDBSCAN_DEFAULT_SETTINGS)
        for key, value in DEFAULT_SETTINGS.items():
            if key not in settings:
                settings[key] = value
        self.settings = settings
        super().add_default_settings()

    @override
    def run_clustering(self) -> None:
        self.init_clusterer()
        self.fit_hdbscan()
        if self.settings["plot_path"] is not None or self.settings["show_plot"] is True:
            self.output_hdbscan()

    def init_clusterer(self):
        self.clusterer = hdbscan.HDBSCAN(
            min_cluster_size=self.settings["min_cluster_size"],
            min_samples=self.settings["min_samples"],
            metric=self.settings["metric"],
            p=self.settings["p"],
            alpha=self.settings["alpha"],
            cluster_selection_epsilon=self.settings["cluster_selection_epsilon"],
            algorithm=self.settings["algorithm"],
            leaf_size=self.settings["leaf_size"],
            approx_min_span_tree=self.settings["approx_min_span_tree"],
            gen_min_span_tree=self.settings["gen_min_span_tree"],
            core_dist_n_jobs=self.settings["core_dist_n_jobs"],
            cluster_selection_method=self.settings["cluster_selection_method"],
            allow_single_cluster=self.settings["allow_single_cluster"],
            prediction_data=self.settings["prediction_data"],
            match_reference_implementation=self.settings["match_reference_implementation"]
        )

    def fit_hdbscan(self) -> None:
        self.labels = self.clusterer.fit_predict(self.points)
        self.clustered = (self.labels >= 0)

    def output_hdbscan(self) -> None:
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(self.settings["width"] / dpi, self.settings["height"] / dpi))
        ax = fig.add_subplot(111)
        ax.set_facecolor(self.settings["background"])

        point_size = 300.0 / np.sqrt(self.points.shape[0])

        ax.scatter(
            self.points[~self.clustered, 0],
            self.points[~self.clustered, 1],
            color=(0.5, 0.5, 0.5),
            s=point_size,
            alpha=0.5
        )
        ax.scatter(
            self.points[self.clustered, 0],
            self.points[self.clustered, 1],
            c=self.labels[self.clustered],
            s=point_size,
            cmap=self.settings["cmap"]
        )
        ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
        if self.settings["plot_path"] is not None:
            ax.figure.savefig(self.settings["plot_path"])
        if self.settings["show_plot"] is True:
            plt.show()