from Bio import SeqIO

from .base_reader import Reader
from ..runner_build import base_logger


class FastaReader(Reader):

    def __init__(self):
        super().__init__()
        self.seq_dict = {}

    @base_logger.prog_log(prog_name="Read FASTA file")
    def read_fasta(self, seq_path: str, seq_type: str = "Haplotype"):
        """
        Read a fasta file and update the dictionary 'seq_dict' with sequence names and sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence,  either "Haplotype" or "Amplicon". Default is "Haplotype".
        """
        with open(seq_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                name, seq = record.description, str(record.seq)
                if seq_type == "Amplicon":
                    name = name.split(";")[0]
                self.seq_dict[name] = seq