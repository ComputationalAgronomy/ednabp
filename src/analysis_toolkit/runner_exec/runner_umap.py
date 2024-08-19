from Bio import SeqIO
import matplotlib.cm
import matplotlib.colors
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import sparse
import subprocess
import tempfile
import umap
from umap.plot import _datashade_points, _themes

from analysis_toolkit.runner_build import (base_runner, utils_sequence)

class UmapRunner(base_runner.SequenceRunner):
    """
    Class for managing UMAP analysis.
    """
    def __init__(self, samplesdata):
        super().__init__(samplesdata)
        self.units2targets = {}
        self.index_list = []

    @base_runner.log_execution("Write UMAP index file", "write_umap.log")
    def run_write(self,
            target_list: list[str],
            target_level: str,
            unit_level: str = "species",
            save_dir: str = ".",
            neighbors: int = 15,
            min_dist: float = 0.1,
            random_state: int = 42,
            calc_dist: bool = True,
            dereplicate_sequence: bool = False,
            sample_id_list: list[str] = [],
        ) -> pd.DataFrame:
        """
        Run the UMAP pipeline and write the index TSV file.
        (UMAP parameters reference: https://umap-learn.readthedocs.io/en/latest/parameters.html)
        Step:
            1. Load sample ID list
            2. Load units2fasta and units2targets dictionaries
            3. Write index FASTA file
            4. Run UMAP
            5. Create index DataFrame (index, sequence ID, unit name)
            6. Update index columns (target, source, UMAP coordinates)
            7. Write index TSV file

        :param target_list: A list of targets to be used (e.g., ["FamilyA", "FamilyB", etc]).
        :param target_level: The taxonomic level of the targets (e.g., family, genus, species).
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

        self._load_units2fasta_units2targets(
            target_list=target_list,
            target_level=target_level,
            unit_level=unit_level,
        )

        self._write_index_fasta(
            aln_index_fasta_path=aln_index_fasta_path,
            dereplicate_sequence=dereplicate_sequence
        )

        self._run_umap(
            fasta_path=aln_index_fasta_path,
            save_dir=save_dir,
            neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state,
            calc_dist=calc_dist
        )

        self._create_index_df()
        self._update_index_columns()
        self.index.to_csv(index_path, sep='\t', index=False)

        self.logger.info(f'Saved index TSV to: {index_path}')

        self.analysis_type = "umap_write"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "target_list": target_list,
                "target_level": target_level,
                "unit_level": unit_level,
                "n_neighbors": neighbors,
                "min_dist": min_dist,
                "random_state": random_state,
                "calc_dist": calc_dist,
                "dereplicate_sequence": dereplicate_sequence,
            }
        )

    @base_runner.log_execution("Plot UMAP results", "plot_umap.log")
    def run_plot(self,
        index_path: str,
        n_unit_threshold: int,
        category: str,
        save_dir: str = ".",
        cmap: str = "rainbow",
        show_legend: bool = True
    ):
        """
        Plot UMAP results based on the specified category (unit, target, or all).

        :param index_path: Path to the pre-created UMAP index file.
        :param n_unit_threshold: Minimum number of sequences for an unit to be included in the analysis. The value should be set equal to UMAP n_neighbors.
        :param category: Column name to group the units by, restricted to 'unit', 'target', or 'all'.
        :param save_dir: The directory to save the PNG files in. Default is the current directory.
        :param cmap: The colormap to use for the plots. Default is "rainbow".
        :param show_legend: Whether to show the legend in the plots. Default is True.
        """
        if category not in ['unit', 'target', 'all']:
            raise ValueError("Invalid category. Must be 'unit', 'target', or 'all'.")

        os.makedirs(save_dir, exist_ok=True)

        self.index = pd.read_csv(index_path, sep='\t')

        self.filtered_index = UmapRunner._filter_index_by_unit_occurrence(self.index, n_unit_threshold)

        self._plot_umap_by_category(
            category=category,
            save_dir=save_dir,
            cmap=cmap,
            show_legend=show_legend
        )

        self.analysis_type = "umap_plot"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "index_path": index_path,
                "category": category,
                "save_dir": save_dir,
                "n_unit_threshold": n_unit_threshold,
                "cmap": cmap,
                "show_legend": show_legend
            }
        )

    def _load_units2fasta_units2targets(self,
            target_list: list[str],
            target_level: str,
            unit_level: str,
            sample_id_list: list[str]
        ) -> tuple[dict[str, str], dict[str, str]]:
        """
        Updates UMAP units to FASTA mapping for a given set of targets.
        It also updates a dictionary linking unit labels to their corresponding target labels.
        """
        for target_name in target_list:
            self._load_units2fasta_dict(
                target_name=target_name,
                target_level=target_level,
                unit_level=unit_level
            )
            self.units2targets.update(dict.fromkeys(list(self.units2fasta.keys()), target_name))

    def _write_index_fasta(self,
            aln_index_fasta_path: str,
            dereplicate_sequence: bool
        ):
        """
        Read a units2fasta dict and output an aligned index FASTA file replacing the sequence IDs with indexes.

        :param fasta_path: Path to the output FASTA file.
        :param index_fasta_path: Path to the output index FASTA file.
        :param aln_index_fasta_path: Path to the output index FASTA file after alignment.
        :param dereplicate_sequence: Whether to dereplicate the sequences.
        """
        temp_dir = tempfile.TemporaryDirectory()
        fasta_path = os.path.join(temp_dir.name, 'umap.fa')
        index_fasta_path = os.path.join(temp_dir.name, "input.fa")

        try:
            utils_sequence.write_fasta(units2fasta_dict=self.units2fasta, save_path=fasta_path, dereplicate=dereplicate_sequence)
 
            with open(fasta_path, 'r') as in_handle, open(index_fasta_path, 'w') as out_handle:
                for i, record in enumerate(SeqIO.parse(in_handle, 'fasta')):
                    index = str(i)
                    unit = record.description.rsplit("-", 1)[0]
                    seq_id = record.description

                    self.index_list.append([index, seq_id, unit])

                    record.id = index
                    record.description = ''
                    record.name = index

                    SeqIO.write(record, out_handle, 'fasta')

            utils_sequence.align_fasta(seq_path=index_fasta_path, aln_path=aln_index_fasta_path)

        finally:
            temp_dir.cleanup()

    def _calc_distmx(self,
            fasta_path: str,
            dist_path: str,
            maxdist: float = 1.0,
            termdist: float = 1.0,
            threads: int = 12
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
        self.logger.info("Calculating distance matrix...")

        cmd = [
            "usearch", "-calc_distmx", fasta_path, "-tabbedout", dist_path,
            "-maxdist", str(maxdist), "-termdist", str(termdist)
        ]
        if threads:
            cmd.extend(["-threads", str(threads)])

        self.logger.info("> Running USEARCH command:", cmd)
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"USEARCH finished. Output distance matrix file saved to: {dist_path}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error occurred during the calculation of the distance matrix: {e.stderr}")

    def _load_sparse_dist_matrix(self, dist_path: str):
        """
        Load a sparse distance matrix from a distance matrix file created by the 'calc_distmx' function.

        :param dist_path: Path to the input distance matrix file.
        :return: Sparse distance matrix as a NumPy array.
        """
        self.matrix = pd.read_csv(dist_path, header=None, sep='\t')
        self.logger.info(f"Loading sparse {max(self.matrix[0])+1} x {max(self.matrix[0])+1} distance matrix from: {dist_path}")

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
        base_map = {'A':[1, 0, 0, 0],
                    'C':[0, 1, 0, 0],
                    'G':[0, 0, 1, 0],
                    'T':[0, 0, 0, 1]}

        one_hot_encoded = []
        for base in sequence:
            one_hot_encoded.extend(base_map.get(base, [0, 0, 0, 0]))

        return one_hot_encoded

    def _load_one_hot_matrix(self, fasta_path: str):
        """
        Read in a aligned FASTA file and output a one-hot encoded matrix.

        :param seq_path: Path to the input FASTA file.
        :return: One-hot encoded matrix as a NumPy array.
        """
        self.logger.info(f"Creating one-hot encoded matrix from: {fasta_path}")

        self.matrix = []
        with open(fasta_path, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                self.matrix.append(UmapRunner._sequence_to_one_hot(record.seq))

    def _fit_umap(self,
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
        self.logger.info(f'Creating UMAP embedding with {neighbors} neighbors...')

        self.reducer = umap.UMAP(
            n_neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state,
            metric="precomputed" if calc_dist else "euclidean"
        )

        self.embedding = self.reducer.fit_transform(self.matrix)

    def _run_umap(self,
            fasta_path: str,
            save_dir: str,
            neighbors: int,
            min_dist: float,
            random_state: int,
            calc_dist: bool
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
        self.index = pd.DataFrame(self.index_list, columns=["index", "seq_id", "unit"])

    def _update_index_target_column(self):
        """
        uses the "unit" column of index pd.DataFrame and uses the 'unit2target' dictionary to map unit labels to target labels.
        If a unit label is not found in the dictionary, it assigns the label "unknown".
        """
        target_labels = []
        for unit in self.index["unit"]:
            if unit in self.units2targets:
                target_labels.append(self.units2targets[unit])
            else:
                target_labels.append("unknown")

        self.index["target"] = target_labels

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
        self.index["umap1"] = self.embedding[:,0]
        self.index["umap2"] = self.embedding[:,1]

    def _update_index_columns(self):
        """
        The steps for updating the index DataFrame with source/target labels and UMAP cordinates.
        """
        if self.units2targets is not None:
            self._update_index_target_column()
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
        symbol_map = ["o", "D", "*", "s", "h", "8", "X", "p"]
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
                    "each sample (size mismatch: {} {})".format(
                        labels.shape[0], points.shape[0]
                    )
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
                    Patch(facecolor=color_key[k], label=k) for k in unique_labels
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
                ax.scatter(points[i, 0], points[i, 1], s=point_size, c=colors[i], marker=m[i], alpha=alpha)
            # ax.scatter(points[:, 0], points[:, 1], s=point_size, c=colors, markers=m, alpha=alpha)


        # Color by values
        elif values is not None:
            if values.shape[0] != points.shape[0]:
                raise ValueError(
                    "Values must have a value for "
                    "each sample (size mismatch: {} {})".format(
                        values.shape[0], points.shape[0]
                    )
                )
            ax.scatter(
                points[:, 0], points[:, 1], s=point_size, c=values, cmap=cmap, alpha=alpha
            )

        # No color (just pick the midpoint of the cmap)
        else:

            color = plt.get_cmap(cmap)(0.5)
            ax.scatter(points[:, 0], points[:, 1], s=point_size, c=color)

        if show_legend and legend_elements is not None:
            ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

        return ax

    def _plot_points(points,
            labels,
            markers,
            cmap,
            show_legend,
            values=None,
            color_key=None,
            background="white",
            width=800,
            height=800
        ):

        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)

        if points.shape[0] <= width * height // 10:
            ax = UmapRunner._matplotlib_points(points, ax, labels, markers, values, color_key, cmap, background, width, height, show_legend)
        else:
            ax = _datashade_points(points, ax, labels, values, color_key, cmap, background, width, height, show_legend)

        ax.set(xticks=[], yticks=[])

        return ax

    def _plot_umap(self,
            png_path: str,
            cmap: str,
            show_legend: bool
        ):
        """
        Plot the UMAP embedding and save the plot as a PNG file.
        """
        points = self.subindex[["umap1", "umap2"]].to_numpy()
        ax = UmapRunner._plot_points(
            points=points,
            labels=self.subindex['unit'],
            markers=self.subindex['source'],
            cmap=cmap,
            show_legend=show_legend
        )
        ax.figure.savefig(png_path, bbox_inches='tight')
        self.logger.info(f"Saved PNG to: {png_path}")
    
        # print('\n> Drawing interactive plot...')
        # p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
        # bokeh.plotting.output_file(html_path)
        # bokeh.plotting.save(p)
        # print(f'Saved plot HTML to: {html_path}')

    def _plot_umap_by_category(self,
            category: str,
            prefix: str,
            save_dir: str,
            cmap: str,
            show_legend: bool
        ) -> None:
        """
        Plot the UMAP embedding and save the plot as a PNG file, grouped by the specified category.
        """
        if category == 'all':
            self.logger.info("Drawing PNG for all units...")
            png_path = os.path.join(save_dir, f"all_umap.png")
            self.subindex = self.filtered_index.copy()
            self._plot_umap(png_path, cmap, show_legend)
            return

        unique_values = np.unique(self.filtered_index[category])
        for value in unique_values:
            self.logger.info(f"Drawing PNG for {value}...")
            png_path = os.path.join(save_dir, f"{value}_umap.png")
            self.subindex = self.filtered_index[self.filtered_index[category] == value]
            self._plot_umap(png_path, cmap, show_legend)