from collections import defaultdict
import os
import tempfile
from typing import Literal

from Bio import AlignIO, SeqIO
import pandas as pd

from ..runner_build import (base_logger, utils_sequence, SeqWriter)
from .seq_hdbscan_clusterer import HDBSCAN_DEFAULT_SETTINGS


class NexusWriter(SeqWriter):

    def __init__(self, sampledata, no_verbose: bool = False):
        super().__init__(sampledata, no_verbose)

    @base_logger.prog_log("Write NEXUS file")
    def write_nexus(self,
            index_path: str,
            species_name: str,
            label_type: Literal["hdbscan", "site"],
            save_dir: str = '.',
            sample_id_list: list[str] | None = None,
            **kwargs
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
        :param kwargs: Additional HDBSCAN settings to override default parameters. Only applied when label_type is 'hdbscan'.
        """
        try:
            self.uniq_seqs2label_freq = {}
            NexusWriter._add_hdbscan_default_settings(self, settings=kwargs)
            self._load_sample_id_list(sample_id_list)
            self._get_files_path(save_dir)
            self._write_spc_seq_files(species_name)
            self._write_nexus_seq_part()
            self._write_nexus_freq_part()
            self.logger.info(f"Saved NEXUS file to: {self.nex_path}")

        finally:
            self._clean_tmp_files()

    def _add_hdbscan_default_settings(self, settings: dict):
        for key, value in settings.items():
            if key in HDBSCAN_DEFAULT_SETTINGS:
                HDBSCAN_DEFAULT_SETTINGS[key] = value
        self.hdbscan_settings = HDBSCAN_DEFAULT_SETTINGS

    def _get_files_path(self, save_dir: str):
        os.makedirs(save_dir, exist_ok=True)
        self.tmp_dir = tempfile.TemporaryDirectory(delete=False)

        self.fa_path = os.path.join(self.tmp_dir.name, f"seq.fa")
        self.uniq_fa_path = os.path.join(self.tmp_dir.name, f"seq_uniq.fa")
        self.uniq_aln_path = os.path.join(self.tmp_dir.name, f"seq_uniq.aln")
        self.nex_path = os.path.join(save_dir, f"haplotype_network.nex")

    def _clean_tmp_files(self):
        self.tmp_dir.cleanup()

    def _write_spc_seq_files(self, species_name):
        self._load_units2fasta_dict(
            taxon_name=species_name,
            taxa_level='species',
            unit_level='species',
        )
        utils_sequence.write_fasta(self.units2fasta, save_path=self.fa_path, dereplicate=False)
        utils_sequence.write_fasta(self.units2fasta, save_path=self.uniq_fa_path, dereplicate=True)
        utils_sequence.align_fasta(seq_file=self.uniq_fa_path, aln_file=self.uniq_aln_path) # TODO(SW): This is a logic issuse. This function should be outside write_nexus_file()

    def _write_nexus_seq_part(self):
        AlignIO.convert(self.uniq_aln_path, "fasta", self.nex_path, "nexus", molecule_type="DNA")

    def _write_nexus_freq_part(self, index_path, species_name, label_type):
        self._get_points_labels(index_path, species_name, label_type)
        self._count_uniq_seq_freq()
        self._assemble_nexus_freq_str()
        with open(self.nex_path, 'a') as nex_handle:
            nex_handle.write(self.freq_string)

    def _get_points_labels(self,
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
            import hdbscan
            points = subindex[["umap1", "umap2"]].to_numpy()
            self.seq_labels = hdbscan.HDBSCAN(
                **self.hdbscan_settings
            ).fit_predict(points)
        elif label_type == 'site':
            self.seq_labels = ["taoyuan" if "taoyuan" in i else "keelung" for i in subindex["seq_id"]]
        else:
            raise ValueError("Label type must be 'hdbscan' or 'site'.")

    def _count_uniq_seq_freq(self) -> str:
        """
        Count the frequency of each label category for each unique sequence in the 'uniq_seqs_path' file based on the 'seqs_path' file.

        :param seqs_path: Path to the input FASTA file containing all sequences.
        :param uniq_seqs_path: Path to the input FASTA file containing only unique sequences.
        """
        with open(self.uniq_fa_path, 'r') as uniq_fa_handle, open(self.fa_path, 'r') as fa_handle:
            uniq_fa_records = list(SeqIO.parse(uniq_fa_handle, 'fasta'))
            fa_records = list(SeqIO.parse(fa_handle, 'fasta'))
            for ufr in uniq_fa_records:
                self.uniq_seqs2label_freq[ufr.name] = defaultdict(int)

                for i, fr in enumerate(fa_records):
                    if str(fr.seq) == str(ufr.seq):
                        self.uniq_seqs2label_freq[ufr.name][self.seq_labels[i]] += 1

    def _assemble_nexus_freq_str(self):
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
