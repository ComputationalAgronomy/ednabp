from abc import ABC, abstractmethod
import numpy as np
import os
import pandas as pd

from analysis_toolkit.runner_build import base_logger
from analysis_toolkit.runner_exec import data_container


def log_execution(prog_name: str, log_file: str):
    def decorator(func):
        def wrapper(*args, **kwargs):
            if "save_dir" in kwargs and os.path.exists(kwargs["save_dir"]):
                fh = base_logger._get_file_handler(os.path.join(kwargs["save_dir"], log_file))
                base_logger.logger.addHandler(fh)
            else:
                pass
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
        self.sample_info = samplesdata.sample_info

    def _load_sample_id_list(self, sample_id_list: str = []):
        if sample_id_list == []:
            self.logger.info(f"No sample ID list specified. Using all {len(self.sample_id_list)} samples.")
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
        self.sample2abundance_dict = {}
        self.abundance_df = pd.DataFrame()

    @log_execution("Write richness or abundance data to CSV", "write_csv.log")
    def run_write(self,
            write_type: str,
            taxa_level: str,
            save_dir: str = ".",
            sample_id_list: list[str] = []
        ):
        """
        Write the richness or abundance data to a CSV file.
        Richness data: number of unique species for each target per sample_id. (e.g. How many different species got in sampleA?)
        Abundance data: sum of read for each taxa per sample_id. (e.g. How many sequences got in sampleA?)

        :param write_type: The type of data to write. Can be "richness" or "abundance".
        :param taxa_level: The name of the level to target (e.g., species, family, etc.).
        :param save_dir: The directory to save the CSV file. Default is the current directory.
        :param sample_id_list: A list of sample IDs to write. Default is an empty list (write all samples).
        """
        if write_type not in ["richness", "abundance"]:
            raise ValueError("Invalid write_type. Must be 'richness' or 'abundance'.")

        os.makedirs(save_dir, exist_ok=True)

        self._load_sample_id_list(sample_id_list)

        for sample_id in self.sample_id_used:
            self._update_abundance_df(sample_id, taxa_level) # columns: target_level, Unit(species), Counts, Sample_id

        self._add_sample_info() # columns: target_level, Unit(species), Counts, Sample_id, Site, Year, Month, Sample

        self._filter_abundance_df()

        if write_type == "richness":
            self.df=self.abundance_df.groupby([taxa_level, "Site", "Year", "Month", "Sample"])["Counts"].nunique().reset_index()
            self.df.to_csv(os.path.join(save_dir, 'Species_richness.csv'), index=False)
        if write_type == "abundance":
            self.df = self.abundance_df.groupby([taxa_level, "Site", "Year", "Month", "Sample"])["Counts"].sum().reset_index()
            self.df.to_csv(os.path.join(save_dir, 'Species_abundance.csv'), index=False)
        
        self.analysis_type = "Write species diversity to csv"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "write_type": write_type,
                "taxa_level": taxa_level,
            }
        )

    def _load_hap_size(self, sample_id: str, hap: str) -> int:
        """
        Get the size of a haplotype for a given sample.

        :param sample_id: The ID of the sample.
        :param hap: The ID of the haplotype.
        :return: The size of the haplotype.
        """
        return int(self.sample_data[sample_id].hap_size[hap])

    def _load_abundance_dict(self, sample_id: str, target_level: str, unit_level: str = "species"):
        """
        Get the abundance of a level for a given sample.

        :param sample_id: The ID of the sample.
        :param level: The name of the level.
        :return: A dictionary mapping level names to abundances.
        """
        self.abundance = {}
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            target_name = level_dict[target_level]
            unit_name = level_dict[unit_level]
            key = (target_name, unit_name)
            if key not in self.abundance:
                self.abundance[key] = 0
            size = self._load_hap_size(sample_id, hap)
            self.abundance[key] += size # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}
    
    def _normalize_abundance(self):
        """
        Normalize the abundance values in the dictionary to percentages.

        :param abundance_dict: A dictionary with unit names and their abundance.
        :returns: A dictionary with unit names and their normalized abundance in percentages.
        """
        total_size = sum(self.abundance.values())
        self.abundance = {key: value/total_size * 100 for key, value in self.abundance.items()}
    
    def _update_sample2abundance_dict(self, sample_id: str):
        self.sample2abundance_dict[sample_id] = self.abundance.copy()
 
    def _abundance_dict2df(self, sample_id, target_level):
        abundance_list = []
        for key, value in self.abundance.items():
            abundance_list.append([key[0], key[1], value])

        self.df = pd.DataFrame(abundance_list, columns=[target_level, "Unit", "Counts"])
        self.df["Sample_id"] = sample_id

    def _update_abundance_df(self, sample_id: str, target_level: str):
        self._load_abundance_dict(sample_id, target_level)
        self._abundance_dict2df(sample_id, target_level)
        self.abundance_df = pd.concat([self.abundance_df, self.df], ignore_index=True)

    def _add_sample_info(self):
        self.abundance_df = pd.merge(left=self.abundance_df, right=self.sample_info, on="Sample_id", how="outer")

    def _filter_abundance_df(self, site_occur_thres: int = 0, sample_occur_thres: int = 0):
        self.abundance_df = self.abundance_df[self.abundance_df['Counts']!=0]

        sp_list = self.abundance_df.groupby(["Unit"])[["Site","Year","Month","Sample"]].nunique()
        sp_list = sp_list[sp_list["Site"]>0] # Filtering by number of sites
        sp_list = sp_list[sp_list["Sample"]>0] # Filtering by number of samples

        survive = np.array([])
        for m in sp_list.index:
            survive = np.append(survive, np.where(self.abundance_df[sp_list.index.name]==m)[0])

        self.abundance_df = self.abundance_df.iloc[survive,:].reset_index(drop=True)
