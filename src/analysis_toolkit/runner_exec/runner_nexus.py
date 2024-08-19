from Bio import AlignIO, SeqIO
import os
import numpy as np
import pandas as pd
import tempfile

from analysis_toolkit.runner_build import (base_runner, utils_sequence)
from analysis_toolkit.runner_exec import runner_hdbscan


class NexusRunner(base_runner.SequenceRunner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)
        self.uniq_seqs2label_freq = {}

    @base_runner.log_execution("Write NEXUS file", "write_nexus.log")
    def run_write(self,
            index_path: str,
            species_name: str,
            label_type: str,
            save_dir: str = '.',
            sample_id_list: list[str] = []
        ) -> None:
        """
        Write a NEXUS file for a given species. The file can be used as input for Popart to plot a haplotype network.
        The NEXUS file contains the aligned sequences, their corresponding labels, and the frequency of each label for each unique sequence.
        'hdbscan' label type uses HDBSCAN clustering results for labeling and calculates the frequency of each cluster label.
        'site' label type uses site information (e.g., 'taoyuan' or 'keelung') for labeling and calculates the frequency of each site label.

        :param index_path: Path to the index file containing the unit information.
        :param species_name: Name of the species for which the NEXUS file will be generated.
        :param label_type: Type of labels to use. Either 'hdbscan' or 'site'.
        :param save_dir: Directory where the NEXUS file will be saved. Default is the current directory.
        :param sample_id_list: List of sample IDs to include in the NEXUS file. The list should be same as that specified by the index file. Default is None (plot all samples).
        """
        try:
            self._load_sample_id_list(sample_id_list)

            self._load_points_labels(
                index_path=index_path,
                species_name=species_name,
                label_type=label_type
            )

            self._load_units2fasta_dict(
                target_name=species_name,
                target_level='species',
                unit_level='species',
            )

            os.makedirs(save_dir, exist_ok=True)
            temp_dir = tempfile.TemporaryDirectory()
            fasta_path = os.path.join(temp_dir.name, f"{species_name}.fa")
            uniq_fasta_path = os.path.join(temp_dir.name, f"{species_name}_uniq.fa")
            aln_fasta_path = os.path.join(temp_dir.name, f"{species_name}_uniq.aln")
            nex_path = os.path.join(save_dir, f"{species_name}.nex")
            self._write_seq_files(fasta_path, uniq_fasta_path, aln_fasta_path, nex_path)

            self._count_uniq_seq_frequency(fasta_path, uniq_fasta_path)
            self._assemble_nex_format_freq_string()
            self._add_freq_string2nex_path(nex_path)

            self.logger.info(f"Saved NEXUS file to: {nex_path}")

            self.analysis_type = "nexus_write"
            self.results_dir = save_dir
            self.parameters.update(
                {
                    "index_path": index_path,
                    "species_name": species_name,
                    "label_type": label_type
                }
            )

        finally:
            temp_dir.cleanup()

    def run_plot(self):
        return super().run_plot()

    def _load_points_labels(self,
            index_path: str,
            species_name: str,
            label_type: str
        ) -> list[str]:
        """
        load points labels from index file.
        """
        index = pd.read_csv(index_path, sep='\t')
        subindex = index[index["unit"] == species_name]
        if label_type == 'hdbscan':
            points = subindex[["umap1", "umap2"]].to_numpy()
            self.seq_labels, _, _, _ = runner_hdbscan.HdbscanRunner._fit_hdbscan(points=points, min_samples=10, min_cluster_size=5) # TODO(SW): FIX THIS
        elif label_type == 'site':
            self.seq_labels = ["taoyuan" if "taoyuan" in i else "keelung" for i in subindex["seq_id"]]
        else:
            raise ValueError("Label type must be 'hdbscan' or 'site'.")

    def _write_seq_files(self, fasta_path, uniq_fasta_path, aln_fasta_path, nex_path):
        utils_sequence.write_fasta(self.units2fasta, save_path=fasta_path, dereplicate=False)
        utils_sequence.write_fasta(self.units2fasta, save_path=uniq_fasta_path, dereplicate=True)
        utils_sequence.align_fasta(seq_file=uniq_fasta_path, aln_file=aln_fasta_path) # TODO(SW): This is a logic issuse. This function should be outside write_nexus_file()
        AlignIO.convert(aln_fasta_path, "fasta", nex_path, "nexus", molecule_type="DNA")

    def _count_uniq_seq_frequency(self, fasta_path, uniq_fasta_path) -> str:
        """
        Count the frequency of each label category for each unique sequence in the 'uniq_seqs_path' file based on the 'seqs_path' file.

        :param seqs_path: Path to the input FASTA file containing all sequences.
        :param uniq_seqs_path: Path to the input FASTA file containing only unique sequences.
        """
        self.uniq_labels = np.unique(self.seq_labels)

        with open(uniq_fasta_path, 'r') as uniq_handle:
            uniq_records = list(SeqIO.parse(uniq_handle, 'fasta'))
            for record in uniq_records:
                self.uniq_seqs2label_freq[record.name] = {uniq_label: 0 for uniq_label in self.uniq_labels}

        with open(fasta_path, 'r') as seq_handle:
            for i, seq_record in enumerate(SeqIO.parse(seq_handle, 'fasta')):
                for uniq_record in uniq_records:
                    if seq_record.seq == uniq_record.seq:
                        self.uniq_seqs2label_freq[uniq_record.name][self.seq_labels[i]] += 1

    def _assemble_nex_format_freq_string(self):
        self.freq_string = (
            f"Begin Traits;\n"
            f"Dimensions NTraits={len(self.uniq_labels)};\n"
            f"Format labels=yes missing=? separator=Comma;\n"
            f"TraitLabels {' '.join(self.uniq_labels)};\n"
            f"Matrix\n"
        )

        for seq_id, label_freq in self.uniq_seqs2label_freq.items():
            freq_values = ",".join(map(str, label_freq.values()))
            self.freq_string += f"{seq_id} {freq_values}\n"

        self.freq_string += ";\nend;\n"
    
    def _add_freq_string2nex_path(self, nex_path):
        with open(nex_path, 'a') as nex_handle:
            nex_handle.write(self.freq_string)

