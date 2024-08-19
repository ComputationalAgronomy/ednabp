from Bio import SeqIO

from analysis_toolkit.read import read_blast_csv
from analysis_toolkit.runner_build import base_logger


class FastaReader(read_blast_csv.Reader):

    def __init__(self):
        super().__init__()
        self.seq_dict = {}

    def parse_fasta(self, seq_path: str, seq_type: str):
        """
        Parse the fasta file and return the dictionary 'seq_dict' containing sequence names and sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence, either "Haplotype" or "Amplicon".
        :return: A dictionary containing sequence names and sequences.
        """
        with open(seq_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                name, seq = record.description, str(record.seq)
                if seq_type == "Amplicon":
                    name = name.split(";")[0]
                self.seq_dict[name] = seq

    def read_fasta(self, seq_path: str, seq_type: str = "Haplotype"):
        """
        Read a fasta file and update the dictionary 'seq_dict' with sequence names and sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence,  either "Haplotype" or "Amplicon". Default is "Haplotype".
        """
        base_logger.logger.info(f"Reading {seq_type} FASTA file: {seq_path}.")

        self.parse_fasta(seq_path, seq_type)

        read_count = len(self.seq_dict)
        base_logger.logger.info(f"COMPLETE: {seq_type} Sequences Read. Total {read_count} reads.")
