import os
import tempfile
from typing import Literal

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib.patches import Patch

from ..runner_build import SeqWriter, base_logger, utils, utils_sequence
from . import seq_hdbscan_clusterer


class UmapRunner(SeqWriter):
    """
    Class for managing UMAP analysis.
    """

    def __init__(self, sampledata, no_verbose=False):
        super().__init__(sampledata, no_verbose)
        self.units2taxa = {}
        self.index_list = []

    @base_logger.prog_log("Write UMAP index file")
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
    ) -> pd.DataFrame:
        """
        Run the UMAP pipeline and write the index TSV file.
        (UMAP parameters reference: https://umap-learn.readthedocs.io/en/latest/parameters.html)
        Step:
            1. Load sample ID list
            2. Load units2fasta and units2taxa dictionaries
            3. Write index FASTA file
            4. Run UMAP
            5. Create index DataFrame (index, sequence ID, unit name)
            6. Update index columns (taxa, source, UMAP coordinates)
            7. Write index TSV file

        :param taxa_list: A list of taxas to be used (e.g., ["FamilyA", "FamilyB", etc]).
        :param taxa_level: The taxonomic level of the taxas (e.g., family, genus, species).
        :param units_level: The taxonomic level of the units. Default is "species"
        :param save_dir: Directory where the output files (FASTA, aligned FASTA, index) will be saved. Default is current directory.
        :param n_neighbors: Number of neighbors to consider for UMAP. Default is 15.
        :param min_dist: Minimum distance parameter for UMAP. Default is 0.1.
        :param random_state: Random seed for reproducibility. Default is 42.
        :param calc_dist: If True, calculates a distance matrix for UMAP. Otherwise, transforms sequences into a one-hot encoded matrix. Default is True.
        :param dereplicate_sequence: If True, use unique sequences as input data for UMAP. Default is False.
        :param sample_id_list: A list of sample IDs to use for UMAP. Default is None (use all samples).
        """
        os.makedirs(save_dir, exist_ok=True)
        index_path = os.path.join(save_dir, "umap_index.tsv")
        aln_index_fasta_path = os.path.join(save_dir, "input.aln")

        self._load_sample_id_list(sample_id_list)

        self._load_units2fasta_units2taxa(
            taxa_list=taxa_list,
            taxa_level=taxa_level,
            unit_level=unit_level,
        )

        self._write_index_fasta(
            aln_index_fasta_path=aln_index_fasta_path,
            dereplicate_sequence=dereplicate_sequence,
        )

        self._run_umap(
            fasta_path=aln_index_fasta_path,
            save_dir=save_dir,
            neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state,
            calc_dist=calc_dist,
        )

        self._create_index_df()
        self._update_index_columns()
        self.index.to_csv(index_path, sep="\t", index=False)
        self.logger.info(f"Index TSV saved_to: {index_path}")

    @base_logger.prog_log("Plot UMAP embedding")
    def plot_umap(
        self,
        index_path: str,
        n_unit_threshold: int,
        category: Literal["taxa", "unit", "all"],
        save_dir: str = ".",
        cmap: str = "rainbow",
        show_legend: bool = True,
    ):
        """
        Plot UMAP results based on the specified category (unit, taxa, or all).

        :param index_path: Path to the pre-created UMAP index file.
        :param n_unit_threshold: Minimum number of sequences for an unit to be included in the analysis. The value should be set equal to UMAP n_neighbors.
        :param category: Column name to group the units by, restricted to 'unit', 'taxa', or 'all'.
        :param save_dir: The directory to save the PNG files in. Default is the current directory.
        :param cmap: The colormap to use for the plots. Default is "rainbow".
        :param show_legend: Whether to show the legend in the plots. Default is True.
        """
        if category not in ["unit", "taxa", "all"]:
            raise ValueError(
                "Invalid category. Must be 'unit', 'taxa', or 'all'."
            )
        self.index = pd.read_csv(index_path, sep="\t")
        self.filtered_index = UmapRunner._filter_index_by_unit_occurrence(
            self.index, n_unit_threshold
        )
        self._plot_umap_by_category(
            category=category,
            save_dir=save_dir,
            cmap=cmap,
            show_legend=show_legend,
        )

    @base_logger.prog_log("Run HDBSCAN clustering for UMAP embedding")
    def hdbscan_umap(
        self,
        index_path: str,
        n_unit_threshold: int,
        category: Literal["taxa", "unit", "all"],
        save_dir=None,
        **settings,
    ) -> pd.DataFrame:
        self.index = pd.read_csv(index_path, sep="\t")
        self.filtered_index = UmapRunner._filter_index_by_unit_occurrence(
            self.index, n_unit_threshold
        )
        return self._hdbscan_umap_by_category(
            category=category, save_dir=save_dir, settings=settings
        )

    def _load_units2fasta_units2taxa(
        self,
        taxa_list: list[str],
        taxa_level: str,
        unit_level: str,
    ) -> tuple[dict[str, str], dict[str, str]]:
        """
        Updates UMAP units to FASTA mapping for a given set of taxas.
        It also updates a dictionary linking unit labels to their corresponding taxa labels.
        """
        for taxon_name in taxa_list:
            self._load_units2fasta_dict(
                taxon_name=taxon_name,
                taxa_level=taxa_level,
                unit_level=unit_level,
            )
            self.units2taxa.update(
                dict.fromkeys(list(self.units2fasta.keys()), taxon_name)
            )

    def _write_index_fasta(
        self, aln_index_fasta_path: str, dereplicate_sequence: bool
    ):
        """
        Read a units2fasta dict and output an aligned index FASTA file replacing the sequence IDs with indexes.

        :param fasta_path: Path to the output FASTA file.
        :param index_fasta_path: Path to the output index FASTA file.
        :param aln_index_fasta_path: Path to the output index FASTA file after alignment.
        :param dereplicate_sequence: Whether to dereplicate the sequences.
        """
        temp_dir = tempfile.TemporaryDirectory()
        fasta_path = os.path.join(temp_dir.name, "umap.fa")
        index_fasta_path = os.path.join(temp_dir.name, "input.fa")

        try:
            utils_sequence.write_fasta(
                units2fasta_dict=self.units2fasta,
                save_path=fasta_path,
                dereplicate=dereplicate_sequence,
            )

            with (
                open(fasta_path) as in_handle,
                open(index_fasta_path, "w") as out_handle,
            ):
                for i, record in enumerate(SeqIO.parse(in_handle, "fasta")):
                    index = str(i)
                    unit = record.description.split("-")[0]
                    seq_id = record.description

                    self.index_list.append([index, seq_id, unit])

                    record.id = index
                    record.description = ""
                    record.name = index

                    SeqIO.write(record, out_handle, "fasta")
            utils_sequence.align_fasta(
                seq_path=index_fasta_path, aln_path=aln_index_fasta_path
            )

        finally:
            temp_dir.cleanup()

    @base_logger.prog_log("Calculate distance matrix")
    def _calc_distmx(
        self,
        fasta_path: str,
        dist_path: str,
        maxdist: float = 1.0,
        termdist: float = 1.0,
        threads: int = 12,
    ):
        """
        Calculate distance matrix using USEARCH.
        (USEARCH command reference: https://drive5.com/usearch/manual/cmd_calc_distmx.html)

        :param seq_path: Path to the input aligned FASTA file.
        :param dist_path: Path to the output distance matrix file.
        :param maxdist: The maximum distance to be written. Default is 1.0.
        :param termdist: The distance threshold for terminating the calculation. Default is 1.0.
        :param threads: Number of threads to use for the calculation. Default is 12.
        """
        cmd = [
            "usearch",
            "-calc_distmx",
            fasta_path,
            "-tabbedout",
            dist_path,
            "-maxdist",
            str(maxdist),
            "-termdist",
            str(termdist),
        ]
        if threads:
            cmd.extend(["-threads", str(threads)])
        utils.run_subprocess("USEARCH", cmd, dist_path)

    @base_logger.prog_log("Load distance matrix")
    def _load_sparse_dist_matrix(self, dist_path: str):
        """
        Load a sparse distance matrix from a distance matrix file created by the 'calc_distmx' function.

        :param dist_path: Path to the input distance matrix file.
        :return: Sparse distance matrix as a NumPy array.
        """
        from scipy import sparse

        self.matrix = pd.read_csv(dist_path, header=None, sep="\t")
        self.logger.info(
            f"Loading sparse {max(self.matrix[0]) + 1} x {max(self.matrix[0]) + 1} distance matrix from: {dist_path}"
        )

        diagonal = self.matrix[0] == self.matrix[1]
        row = np.concatenate([self.matrix[0], self.matrix[1][~diagonal]])
        col = np.concatenate([self.matrix[1], self.matrix[0][~diagonal]])
        data = 1 - np.concatenate([self.matrix[2], self.matrix[2][~diagonal]])

        self.matrix = sparse.csr_matrix((data, (row, col)), dtype=np.float32)

    def _sequence_to_one_hot(sequence: str):
        """
        Convert a sequence with only ATCG bases to a one-hot encoded vector.

        :param sequence: DNA sequence.
        :return: One-hot encoded vector.
        """
        base_map = {
            "A": [1, 0, 0, 0],
            "C": [0, 1, 0, 0],
            "G": [0, 0, 1, 0],
            "T": [0, 0, 0, 1],
        }

        one_hot_encoded = []
        for base in sequence:
            one_hot_encoded.extend(base_map.get(base, [0, 0, 0, 0]))

        return one_hot_encoded

    @base_logger.prog_log("Load one-hot encoded matrix")
    def _load_one_hot_matrix(self, fasta_path: str):
        """
        Read in a aligned FASTA file and output a one-hot encoded matrix.

        :param seq_path: Path to the input FASTA file.
        :return: One-hot encoded matrix as a NumPy array.
        """
        self.logger.info(f"Creating one-hot encoded matrix from: {fasta_path}")

        self.matrix = []
        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.matrix.append(UmapRunner._sequence_to_one_hot(record.seq))

    @base_logger.prog_log("Create UMAP embedding")
    def _fit_umap(
        self,
        neighbors: int,
        min_dist: float,
        random_state: int,
        calc_dist: bool,
    ):
        """
        Fit UMAP and store the UMAP object and the embedding.

        :param neighbors: Number of neighbors for umap.
        :param min_dist: Minimum distance for umap.
        :param random_state: Random state for umap.
        :param precomputed: Whether the elements of the matrix are distances or not.
        """
        import umap

        self.reducer = umap.UMAP(
            n_neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state,
            metric="precomputed" if calc_dist else "euclidean",
        )
        self.embedding = self.reducer.fit_transform(self.matrix)

    def _run_umap(
        self,
        fasta_path: str,
        save_dir: str,
        neighbors: int,
        min_dist: float,
        random_state: int,
        calc_dist: bool,
    ):
        if calc_dist:
            dist_path = os.path.join(save_dir, "distance.txt")
            self._calc_distmx(fasta_path, dist_path)
            self._load_sparse_dist_matrix(dist_path)
        else:
            self._load_one_hot_matrix(fasta_path)

        self._fit_umap(neighbors, min_dist, random_state, calc_dist)

    def _create_index_df(self):
        """
        Creates an index DataFrame with columns for index, sequence ID, unit name, and UMAP coordinates.
        """
        self.index = pd.DataFrame(
            self.index_list, columns=["index", "seq_id", "unit"]
        )

    def _update_index_taxa_column(self):
        """
        uses the "unit" column of index pd.DataFrame and uses the 'unit2taxa' dictionary to map unit labels to taxa labels.
        If a unit label is not found in the dictionary, it assigns the label "unknown".
        """
        taxa_labels = []
        for unit in self.index["unit"]:
            if unit in self.units2taxa:
                taxa_labels.append(self.units2taxa[unit])
            else:
                taxa_labels.append("unknown")

        self.index["taxa"] = taxa_labels

    def _update_index_source_column(self):
        """
        uses the "seq_id" column and identifies if the sequence ID contains the substring recorded in the 'sources' variable and assigns the corresponding label.
        If neither substring is found, it assigns the label "unknown".
        """
        sources = ["taoyuan", "keelung"]

        source_labels = []
        for id in self.index["seq_id"]:
            label = "unknown"
            for source in sources:
                if source in id:
                    label = source
                    break
            source_labels.append(label)

        self.index["source"] = source_labels

    def _updata_index_embedding_columns(self):
        """
        Updates the index DataFrame with UMAP coordinates.
        """
        self.index["umap1"] = self.embedding[:, 0]
        self.index["umap2"] = self.embedding[:, 1]

    def _update_index_columns(self):
        """
        The steps for updating the index DataFrame with source/taxa labels and UMAP cordinates.
        """
        if self.units2taxa is not None:
            self._update_index_taxa_column()
        self._update_index_source_column()
        self._updata_index_embedding_columns()

    @staticmethod
    def _filter_index_by_unit_occurrence(index, n: int = 1) -> pd.DataFrame:
        """
        Filter out units that occur less than n times in the index DataFrame.
        The value of 'n' should be set equal to the value of the UMAP 'neighbor'.

        :param index: Index DataFrame.
        :param n: The threshold of minimum occurrence to keep a unit.
        :return: Filtered index DataFrame.
        """
        if n <= 1:
            return index
        counts = index["unit"].value_counts()
        units_to_remove = counts[counts < n].index
        filtered_index = index[~index["unit"].isin(units_to_remove)]
        return filtered_index

    def _matplotlib_points(
        points,
        ax=None,
        labels=None,
        markers=None,
        values=None,
        color_key=None,
        cmap="rainbow",
        background="white",
        width=800,
        height=800,
        show_legend=True,
        alpha=None,
        symbol_map=("o", "D", "*", "s", "h", "8", "X", "p"),
    ):
        point_size = 300.0 / np.sqrt(points.shape[0])

        legend_elements = None

        # if ax is None:
        #     dpi = plt.rcParams["figure.dpi"]
        #     fig = plt.figure(figsize=(width / dpi, height / dpi))
        #     ax = fig.add_subplot(111)

        ax.set_facecolor(background)

        # Color by labels
        if labels is not None:
            if labels.shape[0] != points.shape[0]:
                raise ValueError(
                    "Labels must have a label for "
                    f"each sample (size mismatch: {labels.shape[0]} {points.shape[0]})"
                )
            if color_key is None:
                unique_labels = np.unique(labels)
                num_labels = unique_labels.shape[0]
                color_key = plt.get_cmap(cmap)(np.linspace(0, 1, num_labels))
                legend_elements = [
                    Patch(facecolor=color_key[i], label=k)
                    for i, k in enumerate(unique_labels)
                ]

            if isinstance(color_key, dict):
                colors = pd.Series(labels).map(color_key)
                unique_labels = np.unique(labels)
                legend_elements = [
                    Patch(facecolor=color_key[k], label=k)
                    for k in unique_labels
                ]
            else:
                unique_labels = np.unique(labels)
                if len(color_key) < unique_labels.shape[0]:
                    raise ValueError(
                        "Color key must have enough colors for the number of labels"
                    )
                new_color_key = {
                    k: matplotlib.colors.to_hex(color_key[i])
                    for i, k in enumerate(unique_labels)
                }
                legend_elements = [
                    Patch(facecolor=color_key[i], label=k)
                    for i, k in enumerate(unique_labels)
                ]
                colors = pd.Series(labels).map(new_color_key)

            if markers is not None:
                m = []
                unique_markers = np.unique(markers)
                if len(unique_markers) > len(symbol_map):
                    raise ValueError(
                        "Too many unique markers for the number of labels, please customize 'symbol_map'."
                    )
                for marker in markers:
                    for i, k in enumerate(unique_markers):
                        if marker == k:
                            m.append(symbol_map[i])
            colors = list(colors)
            for i in range(len(points[:, 0])):
                ax.scatter(
                    points[i, 0],
                    points[i, 1],
                    s=point_size,
                    c=colors[i],
                    marker=m[i],
                    alpha=alpha,
                )
            # ax.scatter(points[:, 0], points[:, 1], s=point_size, c=colors, markers=m, alpha=alpha)

        # Color by values
        elif values is not None:
            if values.shape[0] != points.shape[0]:
                raise ValueError(
                    "Values must have a value for "
                    f"each sample (size mismatch: {values.shape[0]} {points.shape[0]})"
                )
            ax.scatter(
                points[:, 0],
                points[:, 1],
                s=point_size,
                c=values,
                cmap=cmap,
                alpha=alpha,
            )

        # No color (just pick the midpoint of the cmap)
        else:
            color = plt.get_cmap(cmap)(0.5)
            ax.scatter(points[:, 0], points[:, 1], s=point_size, c=color)

        if show_legend and legend_elements is not None:
            ax.legend(
                handles=legend_elements,
                loc="center left",
                bbox_to_anchor=(1, 0.5),
            )

        return ax

    def _plot_points(
        points,
        labels,
        markers,
        cmap,
        show_legend,
        values=None,
        color_key=None,
        background="black",
        width=800,
        height=800,
    ):
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)

        if points.shape[0] <= width * height // 10:
            ax = UmapRunner._matplotlib_points(
                points,
                ax,
                labels,
                markers,
                values,
                color_key,
                cmap,
                background,
                width,
                height,
                show_legend,
            )
        else:
            from umap.plot import _datashade_points

            ax = _datashade_points(
                points,
                ax,
                labels,
                values,
                color_key,
                cmap,
                background,
                width,
                height,
                show_legend,
            )

        ax.set(xticks=[], yticks=[])

        return ax

    def _plot_umap(self, png_path: str, cmap: str, show_legend: bool):
        """
        Plot the UMAP embedding and save the plot as a PNG file.
        """
        points = self.subindex[["umap1", "umap2"]].to_numpy()
        ax = UmapRunner._plot_points(
            points=points,
            labels=self.subindex["unit"],
            markers=self.subindex["source"],
            cmap=cmap,
            show_legend=show_legend,
            background="black",
        )
        ax.figure.savefig(png_path, bbox_inches="tight")
        self.logger.info(f"Saved PNG to: {png_path}")

        # print('\n> Drawing interactive plot...')
        # p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
        # bokeh.plotting.output_file(html_path)
        # bokeh.plotting.save(p)
        # print(f'Saved plot HTML to: {html_path}')

    def _plot_umap_by_category(
        self, category: str, save_dir: str, cmap: str, show_legend: bool
    ) -> None:
        """
        Plot the UMAP embedding and save the plot as a PNG file, grouped by the specified category.
        """
        os.makedirs(save_dir, exist_ok=True)

        if category == "all":
            self.logger.info("Drawing PNG for all units...")
            png_path = os.path.join(save_dir, "all_umap.png")
            self.subindex = self.filtered_index.copy()
            self._plot_umap(png_path, cmap, show_legend)
            return

        unique_values = np.unique(self.filtered_index[category])
        for value in unique_values:
            self.logger.info(f"Drawing PNG for {value}...")
            png_path = os.path.join(save_dir, f"{value}_umap.png")
            self.subindex = self.filtered_index[
                self.filtered_index[category] == value
            ]
            self._plot_umap(png_path, cmap, show_legend)

    def _hdbscan_umap_by_category(
        self, category: str, save_dir, settings: dict
    ) -> pd.DataFrame:
        df = pd.DataFrame(
            columns=[
                "class",
                "actual_num",
                "cluster_num",
                "cluster_perc",
                "silhouette_avg",
                "ari",
            ]
        )
        if category == "all":
            self.logger.info("HDBSCAN UMAP embeddings for all units...")
            points = self.filtered_index[["umap1", "umap2"]].to_numpy()
            true_labels = self.filtered_index["unit"]
            df.loc[len(df)] = ["all"] + list(
                seq_hdbscan_clusterer.HdbClusterer().run(
                    points=points, true_labels=true_labels, **settings
                )
            )
            return df

        unique_values = np.unique(self.filtered_index[category])
        for value in unique_values:
            self.logger.info(f"HDBSCAN UMAP embeddings for {value}...")
            self.subindex = self.filtered_index[
                self.filtered_index[category] == value
            ]
            points = self.subindex[["umap1", "umap2"]].to_numpy()
            true_labels = self.subindex["unit"]
            if save_dir:
                plot_path = os.path.join(save_dir, f"hdbscan_{value}.png")
            df.loc[len(df)] = [value] + list(
                seq_hdbscan_clusterer.HdbClusterer().run(
                    points=points,
                    true_labels=true_labels,
                    plot_path=plot_path,
                    **settings,
                )
            )
        return df
