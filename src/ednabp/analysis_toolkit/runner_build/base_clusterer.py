from abc import abstractmethod, ABC
import logging
import os

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score, adjusted_rand_score

class Clusterer(ABC):

    def __init__(
            self,
            index: pd.DataFrame,
        ):
        self.index = index
        self.get_embeddings()
        self.run_clustering()
        self.calc_metrics()
        self.output_metrics()

    def add_default_settings(self) -> None:
        defaults = {
            "show_log": True,
            "metrics_log_path": None,
            "log_overwrite": True,
            "cmap": "Spectral",
            "background": "white",
            "width": 800,
            "height": 800,
        }
        for key, value in defaults.items():
            if key not in self.settings:
                self.settings[key] = value

    # def load_index(self, index_path: str) -> None:
    #     self.index = pd.read_csv(index_path, sep='\t')

    def get_embeddings(self) -> None:
        self.points = self.index[["umap1", "umap2"]].to_numpy()

    @abstractmethod
    def run_clustering(self) -> None:
        pass

    def calc_metrics(self) -> None:
        true_labels = self.index["unit"]
        self.actual_num = len(true_labels.unique())
        self.cluster_num = max(self.labels) + 1

        noise_counts = sum(1 for i in self.labels if i < 0)
        self.cluster_perc = (1 - noise_counts / len(self.labels)) * 100

        if len(np.unique(self.labels)) == 1:
            self.silhouette_avg = 0
            self.ari = 0
        else:
            self.silhouette_avg = silhouette_score(self.points, self.labels)
            self.ari = adjusted_rand_score(true_labels, self.labels)

    def output_metrics(self) -> None:
        FORMAT = "%(message)s"
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(FORMAT) #, TIME_FORMAT)

        if self.settings["show_log"] is not True and self.settings["metrics_log_path"] is None:
            raise ValueError("At least one of `show_log` or `metrics_log_path` must be True")

        if self.settings["show_log"] is True:
            sh = logging.StreamHandler()
            sh.setLevel(logging.INFO)
            sh.setFormatter(formatter)
            logger.addHandler(sh)

        if self.settings["metrics_log_path"] is not None:
            if os.path.exists(self.settings["metrics_log_path"]) and self.settings["log_overwrite"]:
                os.remove(self.settings["metrics_log_path"])
            fh = logging.FileHandler(filename=self.settings["metrics_log_path"])
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