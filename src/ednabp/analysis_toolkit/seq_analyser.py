from typing import TYPE_CHECKING, Literal, Union

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

    from .mito_data import MitoData
    from .sample_data import SampleData


class SeqAnalyser:
    def __init__(
        self,
        sample_data: Union["SampleData", "MitoData", None] = None,
        verbose: bool = True,
    ):
        self.sample_data = sample_data
        self.verbose = verbose

    def _import_mltreewriter(self):
        from .runner_exec.seq_mltree_writer import MLTreeWriter

        self.mltree_writer = MLTreeWriter(self.sample_data, self.verbose)

    def write_mltree(
        self,
        save_dir: str,
        taxa_list: list[str],
        taxa_level: str,
        unit_level: str = "species",
        save_prefix: str = "ml_tree",
        model: str = None,
        bootstrap: int = None,
        threads: int = None,
        dereplicate_sequence: bool = True,
        n_unit_threshold: int = 1,
        sample_id_list: list[str] | None = None,
    ) -> None:
        if not hasattr(self, "mltree_writer"):
            self._import_mltreewriter()
        self.mltree_writer.write_mltree(
            save_dir,
            taxa_list,
            taxa_level,
            unit_level,
            save_prefix,
            model,
            bootstrap,
            threads,
            dereplicate_sequence,
            n_unit_threshold,
            sample_id_list,
        )

    def _import_nexuswriter(self):
        from .runner_exec.seq_nexus_writer import NexusWriter

        self.nexus_writer = NexusWriter(self.sample_data, self.verbose)

    def write_nexus(
        self,
        index_path: str,
        species_name: str,
        label_type: Literal["hdbscan", "site"],
        save_dir: str = ".",
        sample_id_list: list[str] | None = None,
        **kwargs,
    ) -> None:
        if not hasattr(self, "nexus_writer"):
            self._import_nexuswriter()
        self.nexus_writer.write_nexus(
            index_path,
            species_name,
            label_type,
            save_dir,
            sample_id_list,
            **kwargs,
        )

    def _import_umaprunner(self):
        from .runner_exec.seq_umap_runner import UmapRunner

        self.umap_runner = UmapRunner(self.sample_data, self.verbose)

    def write_umap_index(
        self,
        taxa_list: list[str],
        taxa_level: str,
        unit_level: str = "species",
        save_dir: str = ".",
        neighbors: int = 15,
        min_dist: float = 0.1,
        random_state: int = 42,
        calc_dist: bool = True,
        dereplicate_sequence: bool = False,
        sample_id_list: list[str] | None = None,
    ) -> "pd.DataFrame":
        if not hasattr(self, "umap_runner"):
            self._import_umaprunner()
        return self.umap_runner.write_umap_index(
            taxa_list,
            taxa_level,
            unit_level,
            save_dir,
            neighbors,
            min_dist,
            random_state,
            calc_dist,
            dereplicate_sequence,
            sample_id_list,
        )

    def plot_umap(
        self,
        index_path: str,
        n_unit_threshold: int,
        category: Literal["taxa", "unit", "all"],
        save_dir: str = ".",
        cmap: str = "rainbow",
        show_legend: bool = True,
    ) -> None:
        if not hasattr(self, "umap_runner"):
            self._import_umaprunner()
        self.umap_runner.plot_umap(
            index_path, n_unit_threshold, category, save_dir, cmap, show_legend
        )

    def hdbscan_umap(
        self,
        index_path: str,
        n_unit_threshold: int,
        category: Literal["taxa", "unit", "all"],
        save_dir=None,
        **settings,
    ) -> "pd.DataFrame":
        if not hasattr(self, "umap_runner"):
            self._import_umaprunner()
        self.umap_runner.plot_umap(
            index_path, n_unit_threshold, category, save_dir, **settings
        )
