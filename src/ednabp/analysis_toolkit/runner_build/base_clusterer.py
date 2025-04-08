import logging
import os
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, silhouette_score


class Clusterer(ABC):
    def run(
        self,
        index_path: str | None = None,
        points: np.ndarray | None = None,
        true_labels: np.ndarray | None = None,
    ) -> tuple[int, int, float, float, float]:
        self._get_data(index_path, points, true_labels)
        self._run_clustering()
        self._calc_metrics()
        return (
            self.actual_num,
            self.cluster_num,
            self.cluster_perc,
            self.silhouette_avg,
            self.ari,
        )

    def _add_default_settings(self) -> None:
        # default output settings
        defaults = {
            "print_metrics_log": True,
            "metrics_log_path": None,
            "show_plot": True,
            "plot_path": None,
            "cmap": "Spectral",
            "background": "white",
            "width": 800,
            "height": 800,
        }
        for key, value in defaults.items():
            if key not in self.settings:
                self.settings[key] = value

    def _get_data(self, index_path, points, true_labels) -> None:
        if index_path is not None:
            self._load_index(index_path)
            self._get_embeddings()
            self._get_true_labels()
        elif points is not None and true_labels is not None:
            self.points = points
            self.true_labels = true_labels
        else:
            raise ValueError(
                "Either index_path or (points and true_labels) must be provided"
            )

    def _load_index(self, index_path: str) -> None:
        assert os.path.exists(index_path), (
            f"Index file not found: {index_path}"
        )
        self.index = pd.read_csv(index_path, sep="\t")

    def _get_embeddings(self) -> None:
        self.points = self.index[["umap1", "umap2"]].to_numpy()

    def _get_true_labels(self) -> None:
        self.true_labels = self.index["unit"]

    @abstractmethod
    def _run_clustering(self) -> None:
        pass

    def _calc_metrics(self) -> None:
        self.actual_num = len(np.unique(self.true_labels))
        self.cluster_num = max(self.cluster_labels) + 1

        noise_counts = sum(1 for i in self.cluster_labels if i < 0)
        self.cluster_perc = (1 - noise_counts / len(self.cluster_labels)) * 100

        if len(np.unique(self.cluster_labels)) == 1:
            self.silhouette_avg = 0
            self.ari = 0
        else:
            self.silhouette_avg = silhouette_score(
                self.points, self.cluster_labels
            )
            self.ari = adjusted_rand_score(
                self.true_labels, self.cluster_labels
            )

        if (
            self.settings["print_metrics_log"] is not True
            and self.settings["metrics_log_path"] is None
        ):
            return

        FORMAT = "%(message)s"
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(FORMAT)  # , TIME_FORMAT)

        if self.settings["print_metrics_log"] is True:
            sh = logging.StreamHandler()
            sh.setLevel(logging.INFO)
            sh.setFormatter(formatter)

        if self.settings["metrics_log_path"] is not None:
            if os.path.exists(self.settings["metrics_log_path"]):
                os.remove(self.settings["metrics_log_path"])
            fh = logging.FileHandler(
                filename=self.settings["metrics_log_path"]
            )
            fh.setLevel(logging.INFO)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

        logger.info(f"Actual number of units: {self.actual_num}")
        logger.info(f"Number of clusters: {self.cluster_num}")
        logger.info(f"Cluster percentage: {self.cluster_perc:.2f}%")
        logger.info(f"Silhouette score: {self.silhouette_avg:.4f}")
        logger.info(f"Adjusted Rand Index: {self.ari:.4f}")

        while logger.hasHandlers():
            logger.removeHandler(logger.handlers[0])
