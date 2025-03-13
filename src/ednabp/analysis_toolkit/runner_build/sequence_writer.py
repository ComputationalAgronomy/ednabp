from abc import ABC
from collections import defaultdict

from . import base_writer


class SeqWriter(base_writer.Writer, ABC):
    """
    An abstract class for running sequence related analysis
    """

    def __init__(self, sampledata, no_verbose):
        super().__init__(sampledata, no_verbose)
        self.units2fasta = defaultdict(str)

    def _filter_sequence(self, n_unit_threshold):
        for unit, fasta in self.units2fasta.copy().items():
            seq_num = fasta.count(">")
            if seq_num < n_unit_threshold:
                del self.units2fasta[unit]

    def _load_units2fasta_dict(self, **kargs):
        for sample_id in self.sample_id_used:
            for hap, level_dict in self.sample_data[
                sample_id
            ].hap2level.items():
                if kargs["taxon_name"] != level_dict[kargs["taxa_level"]]:
                    continue
                unit_name = level_dict[kargs["unit_level"]]
                title = f"{unit_name}-{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]
                self.units2fasta[unit_name] += f">{title}\n{seq}\n"

        if "n_unit_threshold" in kargs and kargs["n_unit_threshold"] > 0:
            self._filter_sequence(kargs["n_unit_threshold"])
