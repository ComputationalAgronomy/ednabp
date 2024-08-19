from abc import abstractmethod
import csv
import os
import pickle
from datetime import date

from analysis_toolkit.read import read_blast_csv
from analysis_toolkit.read import read_denoise_report
from analysis_toolkit.read import read_fasta
from analysis_toolkit.runner_build import base_logger


class OneSampleData():
    """
    Container for handling and processing data from various bioinformatics files.
    This class initializes by reading data from given FASTA, denoise report, and BLAST table files,
    and stores the parsed information for further use.

    :param uniq_fasta_path: Path to the unique amplicon FASTA file.
    :param zotu_fasta_path: Path to the ZOTU haplotype FASTA file.
    :param denoise_report_path: Path to the denoise report file.
    :param blast_table_path: Path to the BLAST table file.

    :attribute amp_seq: A dictionary containing amplicon sequences from the unique FASTA file.
    :attribute hap_seq: A dictionary containing haplotype sequences from the ZOTU FASTA file.
    :attribute amp_size: A dictionary containing amplicon sizes from the denoise report.
    :attribute hap2amp: A dictionary mapping haplotypes to amplicons from the denoise report.
    :attribute hap_size: A dictionary containing haplotype sizes from the denoise report.
    :attribute hap2level: A dictionary mapping haplotypes to taxonomic levels from the BLAST table.

    """
    def __init__(self,
            uniq_fasta: str,
            denoise_fasta: str,
            denoise_report: str,
            blast_table: str
        ):
        ufr = read_fasta.FastaReader()
        ufr.read_fasta(seq_path=uniq_fasta, seq_type="Amplicon")
        self.amp_seq = ufr.seq_dict

        dfr = read_fasta.FastaReader()
        dfr.read_fasta(seq_path=denoise_fasta, seq_type="Haplotype")
        self.hap_seq = dfr.seq_dict

        drr = read_denoise_report.DenoiseReportReader()
        drr.read_denoise_report(denoise_report=denoise_report)
        self.amp_size = drr.amp_size
        self.hap2amp = drr.hap2amp
        self.hap_size = drr.hap_size

        br = read_blast_csv.BlastReader()
        br.read_blast_table(blast_table=blast_table)
        self.hap2level = br.hap2level


class SampleData():
    """
    A class for managing sample data storage.
    It provides methods for importing, saving, loading.

    :attribute DATA_FILE_INFO: A dictionary mapping import data types to tuples containing the corresponding child directory and file suffix.
    :attribute sample_data: A dictionary to store sample data, the key is sample ID, and the value is a OneSampleData instance.
    :attribute sample_info: A dictionary to store sample information, the key is sample ID, and the value is a dictionary containing sample metadata.
    :attribute sample_id_list: A list to store sample IDs.
    :attribute verbose: A boolean flag to control logging verbosity. Default is True.
    """
    # keys are data type, values are tuples of (child directory, file suffix).
    DATA_FILE_INFO = {
        "uniq_fasta": ("dereplicate", "_uniq.fasta"),
        "denoise_fasta": ("denoise", "_denoise.fasta"),
        "denoise_report": ("denoise", "_denoise_report.txt"),
        "blast_table": ("blast", "_blast.csv")
    }

    def __init__(self,
        verbose = True,
        logger = base_logger.logger
        ):
        self.sample_data = {}
        self.sample_id_list = []
        self.sample_info = {}
        self.verbose = verbose
        self.logger = logger

        if not self.verbose:
            self.logger.setLevel("WARNING")

    def _check_dir(self) -> None:
        """
        Check if all necessary child directories exist within the specified parent directory.
        """
        self.logger.info("Start checking whether directories exist...")
        for child_dir_name, _ in SampleData.DATA_FILE_INFO.values():
            child_dir = os.path.join(self.import_dir, child_dir_name)

            if not os.path.isdir(child_dir):
                raise FileNotFoundError(f"Directory does not exist: {child_dir}.")

        self.logger.info("All directories exist.")

    def _get_sample_id_list(self):
        """
        Retrieve a list of sample IDs from the specified parent directory.

        :param parent_dir: Path to the parent directory containing sample data.
        :return: A list of sample IDs.
        """
        child_dir, suffix = tuple(SampleData.DATA_FILE_INFO.values())[0]

        child_dir_path = os.path.join(self.import_dir, child_dir)

        self.logger.info(f"Searching sample IDs: prefix with the suffix '{suffix}' in the directory: {child_dir_path}.")

        file_list = os.listdir(child_dir_path)
        sample_id_list = [file.replace(suffix, '') for file in file_list if file.endswith(suffix)]
        self.sample_id_list.extend(sample_id_list)

    def _get_file_paths(self, sample_id: str) -> dict[str, str]:
        """
        Retrieve the file paths for all data types for a given sample ID.

        :param sample_id: The sample ID to retrieve file paths for.
        """
        self.file_paths = {}
        for file_key, (child_dir, suffix) in SampleData.DATA_FILE_INFO.items():
            file_path = os.path.join(self.import_dir, child_dir, f"{sample_id}{suffix}")

            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"File does not exist: {file_path}.")

            self.file_paths[file_key] = file_path

    def _read_sample_info(self, sample_info_path: str) -> None:
        """
        Read sample information from a specified file path.
        """
        prog_name = f"Read sample information from: {sample_info_path}."
        self.logger.info(f"Program: {prog_name}")
        with open(sample_info_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                key = row.pop(list(row.keys())[0])
                if key in self.sample_id_list:
                    self.sample_info[key] = dict(row)
        self.logger.info(f"COMPLETE: {prog_name}")
        # {'sample1': {'site': 'siteA', 'date': '05/20/2020', 'sample': '1'},
        #  'sample2': {'site': 'siteA', 'date': '05/20/2020', 'sample': '2'}}

    def _save_instance(self) -> None:
        """
        Save an instance of SamplesContainer to a specified path.
        """
        with open(self.save_instance_path, 'wb') as f:
            pickle.dump(self, f)

    def import_data(self, import_dir: str, sample_id_list: list[str] = [], sample_info_path: str = None) -> None:
        """
        Import sample data from the specified parent directory.
        The parent directory is expected to contain four data types recorded in 'DATA_FILE_INFO'.
        Data should be organized in child directories with specific suffixes as defined in 'DATA_FILE_INFO'.

        :param import_dir: Path to the parent directory containing the sample data.
        :param sample_id_list: List of sample IDs to import. If not provided, all available sample IDs will be imported. The sample IDs are extracted from the file names using the provided suffix.
        :param sample_info_path: Path to the sample information CSV file. If provided, sample information will be loaded from this file. Default is None.
        """
        self.import_dir = import_dir

        prog_name = f"Read samples from: {self.import_dir}."
        self.logger.info(f"Program: {prog_name}")
        self._check_dir()

        if sample_id_list == []:
            self.logger.info("No sample id list provided.")
            self._get_sample_id_list()
        else:
            self.logger.info("Specified sample id list.")
            for sample_id in sample_id_list:
                if sample_id not in self.sample_id_list:
                    self.sample_id_list.append(sample_id)
                else:
                    self.logger.warning(f"WARNING: Duplicate sample ID: {sample_id}.")

        for sample_id in self.sample_id_list:
            self._get_file_paths(sample_id)
            self.sample_data[sample_id] = OneSampleData(**self.file_paths)

        self.logger.info(f"COMPLETE: {prog_name}")

        if sample_info_path:
            self._read_sample_info(sample_info_path)

    def save_data(
            self, save_instance_dir: str,
            save_prefix: str = f"eDNA_samples_{date.today()}",
            overwrite: bool = False
        ):
        """
        Save the current sample data to a specified directory.

        :param save_instance_dir: If provided, save the SamplesContainer instance to the specified directory. Defaults to None.
        :param save_prefix: The prefix for the save .pkl file. Defaults is 'eDNA_samples_<current_date>'.
        :param overwrite: If True, overwrite the existing file. Defaults to False.
        """
        os.makedirs(save_instance_dir, exist_ok=True)

        self.save_instance_path = os.path.join(save_instance_dir, f"{save_prefix}.pkl")
        prog_name = f"Save sample data to: {self.save_instance_path}."
        self.logger.info(f"Program: {prog_name}")

        if os.path.exists(self.save_instance_path) and not overwrite:
            self.logger.warning(f"WARNING: File already exists: {self.save_instance_path}. Data didn't saved.")
            return
        else:
            pass

        self._save_instance()
        self.logger.info(f"COMPLETE: {prog_name}")

    def load_data(self, load_instance_path) -> None:
        """
        Load sample data from a specified path.
    
        :param load_instance_path: If provided, load a pre-existing SamplesContainer instance from the path. Defaults is None.
        """
        prog_name = f"Load sample data from: {load_instance_path}."
        self.logger.info(f"Program: {prog_name}")
        
        with open(load_instance_path,'rb') as file:
            self.__dict__ = pickle.load(file).__dict__
        self.load_instance_path = load_instance_path
        self.logger.info(f"COMPLETE: {prog_name}")

    def merge_data(self, *object_names: object) -> None:
        """
        Merge sample data from another SamplesContainer instance into the current instance.

        :param samplesdata: The SamplesContainer instance to merge data from.
        """
        for object in object_names:
            prog_name = f"Merge object {object} into the current instance."
            self.logger.info(f"Program: {prog_name}")

            if not isinstance(object, SampleData):
                raise TypeError("FAIL: Input object must be an instance of SamplesData.")

            for sample_id in object.sample_id_list:
                if sample_id in self.sample_data.keys():
                    self.logger.warning(f"WARNING: Sample ID {sample_id} already exists in the current instance. Skipping merge.")
                else:
                    self.sample_data[sample_id] = object.sample_data[sample_id]
                    self.sample_info[sample_id] = object.sample_info[sample_id]

            self.sample_id_list.extend(object.sample_id_list)
            self.sample_id_list = list(set(self.sample_id_list))

            self.logger.info(f"COMPLETE: {prog_name}")

