from abc import ABC, abstractmethod
import os

from analysis_toolkit.runner_build import base_logger
from analysis_toolkit.runner_exec import data_container


def log_execution(prog_name: str, log_file: str):
    def decorator(func):
        def wrapper(*args, **kwargs):
            if kwargs["save_dir"]:
                base_logger._add_file_handler(os.path.join(kwargs["save_dir"], log_file))
            base_logger.logger.info(f"Program: {prog_name}")
            result = func(*args, **kwargs)
            base_logger.logger.info(f"COMPLETE: {prog_name}")
            return result
        return wrapper
    return decorator


class Runner(ABC):
    def __init__(self, sampledata: data_container.SampleData, logger=base_logger.logger):
        self.logger = logger
        self.analysis_type = None
        self.sample_id_used = None
        self.parameters = {}
        self.results_dir = None
        self._import_data(sampledata)

    def _import_data(self, samplesdata):
        self.sample_data = samplesdata.sample_data
        self.sample_id_list = samplesdata.sample_id_list

    def _load_sample_id_list(self, sample_id_list: str = []):
        if sample_id_list == []:
            self.logger.info(f"No sample ID list specified. Using all {len(sample_id_list)} samples.")
            self.sample_id_used = self.sample_id_list
        else:
            self.logger.info(f"Specified {len(sample_id_list)} samples.")
            for sample_id in sample_id_list:
                if sample_id not in self.sample_id_list:
                    raise ValueError(f"Specified invalid sample ID: {sample_id}.")
            self.sample_id_used = sample_id_list

    def _add_file_handler(self, log_path):
        fh = base_logger._get_file_handler(log_path)
        self.logger.addHandler(fh)

    @abstractmethod
    def run_write(self):
        pass

    @abstractmethod
    def run_plot(self):
        pass


class SequenceRunner(Runner):
    def __init__(self, sampledata: data_container.SampleData):
        super().__init__(sampledata)
        self.units2fasta = {}

    def _filter_sequence(self, n_unit_threshold):
        for unit, fasta in self.units2fasta.copy().items():
            seq_num = fasta.count('>')
            if seq_num < n_unit_threshold:
                del self.units2fasta[unit]

    def _load_units2fasta_dict(self,
            target_name: str,
            target_level: str,
            unit_level: str,
            n_unit_threshold: int = -1 # TODO(SW): OR *args
        ):
        for sample_id in self.sample_id_used:
            for hap, level_dict in self.sample_data[sample_id].hap2level.items():
                if target_name not in level_dict[target_level]:
                    continue

                unit_name = level_dict[unit_level]
                title = f"{unit_name}-{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in self.units2fasta:
                    self.units2fasta[unit_name] = ""
                self.units2fasta[unit_name] += f'>{title}\n{seq}\n'

        if n_unit_threshold > 1:
            self._filter_sequence(n_unit_threshold)


class AbundanceRunner(Runner):
    def __init__(self, sampledata: data_container.SampleData):
        super().__init__(sampledata)
        self.samples2abundance = {}

    def _load_hap_size(self, sample_id: str, hap: str) -> int:
        """
        Get the size of a haplotype for a given sample.

        :param sample_id: The ID of the sample.
        :param hap: The ID of the haplotype.
        :return: The size of the haplotype.
        """
        return int(self.sample_data[sample_id].hap_size[hap])

    def _load_units2abundance_dict(self, sample_id: str, unit_level: str):
        """
        Get the abundance of a level for a given sample.

        :param sample_id: The ID of the sample.
        :param level: The name of the level.
        :return: A dictionary mapping level names to abundances.
        """
        self.units2abundance = {}
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            level_name = level_dict[unit_level]
            if level_name not in self.units2abundance:
                self.units2abundance[level_name] = 0
            size = self._load_hap_size(sample_id, hap)
            self.units2abundance[level_name] += size # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}
    
    def _normalize_abundance(self):
        """
        Normalize the abundance values in the dictionary to percentages.

        :param abundance_dict: A dictionary with unit names and their abundance.
        :returns: A dictionary with unit names and their normalized abundance in percentages.
        """
        total_size = sum(self.units2abundance.values())
        self.units2abundance = {key: value/total_size * 100 for key, value in self.units2abundance.items()}
    
    def _update_samples2abundance_dict(self, sample_id: str):
        self.samples2abundance[sample_id] = self.units2abundance.copy()