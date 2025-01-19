from datetime import date
import os
import pickle
import re

import pandas as pd

from ..read import taxa_table_reader
from ..read import denoise_report_reader
from ..read import fasta_reader
from ..runner_build import base_logger


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
        ufr = fasta_reader.FastaReader()
        ufr.read_fasta(seq_path=uniq_fasta, seq_type="Amplicon")
        self.amp_seq = ufr.seq_dict

        dfr = fasta_reader.FastaReader()
        dfr.read_fasta(seq_path=denoise_fasta, seq_type="Haplotype")
        self.hap_seq = dfr.seq_dict

        drr = denoise_report_reader.DenoiseReportReader()
        drr.read_denoise_report(denoise_report=denoise_report)
        self.amp_size = drr.amp_size
        self.hap2amp = drr.hap2amp
        self.hap_size = drr.hap_size

        br = taxa_table_reader.TaxaTableReader()
        br.read_taxa_table(blast_table=blast_table)
        self.hap2level = br.hap2level


class SampleData():
    """
    A class for managing sample data storage.
    It provides methods for importing, pickling, unpickling, merging.

    :attribute DATA_FILE_INFO: A dictionary mapping import data types to tuples containing the corresponding child directory and file suffix.
    :attribute sample_data: A dictionary to store sample data, the key is sample ID, and the value is a OneSampleData instance.
    :attribute sample_metadata A dictionary to store sample information, the key is sample ID, and the value is a dictionary containing sample metadata.
    :attribute sample_id_list: A list to store sample IDs.
    :attribute no_verbose: A boolean flag to control logging verbosity. Default is True.
    """

    def __init__(self, no_verbose = False):
        self.sample_data = {}
        self.sample_id_list = []
        self.no_verbose = no_verbose
        self.logger = base_logger.logger

        if self.no_verbose:
            self.logger.setLevel("WARNING")

    def import_data(self,
            dereplicate_dir: str,
            denoise_dir: str,
            assigntaxa_dir: str,
            sample_id_list: list[str] | None = None,
            sample_metadata_path: str | None = None,
            **suffixes
        ):
        """
        Import sample data by specifying the import directories, import file suffixes, and sample metadata table.

        :param dereplicate_dir: Path to the directory containing unique amplicon FASTA files.
        :param denoise_dir: Path to the directory containing denoised haplotype FASTA files and denoise report TXT files.
        :param assigntaxa_dir: Path to the directory containing taxa table CSV files.
        :param sample_id_list: List of sample IDs to import. If not provided, all available sample IDs will be imported. The sample IDs are extracted from the file names using the provided suffix.
        :param sample_metadata_path: Path to the sample metadata CSV file. If provided, sample information will be loaded from this file. Default is None.
        :param suffixes: Additional optional arguments to specify the suffixes format of import files. These include:
            - dereplicate_suffix: Suffix for unique amplicon FASTA files. Default: "_uniq.fasta".
            - denoise_suffix: Suffix for denoised haplotype FASTA files. Default: "_zotus.fasta".
              p.s. Suffix of denoise table would be `re.sub(r"\.\w+", "_report.txt", DENOISE_SUFFIX)`.
            - assigntaxa_suffix: Suffix for taxa CSV table files. Default: "_taxa.csv".
        """
        self._parse_import_info(dereplicate_dir, denoise_dir, assigntaxa_dir, suffixes)

        self._read_sample_data(sample_id_list)

        if sample_metadata_path is not None:
            self._read_sample_metadata(sample_metadata_path)

    def _parse_import_info(self,
            dereplicate_dir: str,
            denoise_dir: str,
            assigntaxa_dir: str,
            suffixes: dict
        ):
        self.suffixes = suffixes
        self._add_default_suffixes()
        self.import_info = [(dereplicate_dir, self.suffixes["dereplicate_suffix"]),
                            (denoise_dir, self.suffixes["denoise_suffix"]),
                            (denoise_dir, self.suffixes["denoise_report_suffix"]),
                            (assigntaxa_dir, self.suffixes["assigntaxa_suffix"])]
        self._check_import_dir()

    def _add_default_suffixes(self):
        DEFAULT_SUFFIXES = {
            "dereplicate_suffix": "_uniqs.fasta",
            "denoise_suffix": "_zotus.fasta",
            "assigntaxa_suffix": "_taxa.csv"
        }
        for key, value in DEFAULT_SUFFIXES.items():
            if key not in self.suffixes:
                self.suffixes[key] = value
        self.suffixes["denoise_report_suffix"] = re.sub(r"\.\w+", "_report.txt", self.suffixes["denoise_suffix"])

    def _check_import_dir(self):
        for in_dir, _ in self.import_info:
            if not os.path.isdir(in_dir):
                raise FileNotFoundError(f"Directory does not exist: {in_dir}.")

    @base_logger.prog_log(prog_name="Import sample sequence data")
    def _read_sample_data(self, sample_id_list: list[str] | None):
        base_logger.print_space(self.logger)
        if sample_id_list is None:
            self.logger.info("No sample id list provided.")
            self._add_unspecified_sample_id_list()
        else:
            self.logger.info("Specified sample id list.")
            self._add_specified_sample_id_list(sample_id_list)

        for sample_id in self.sample_id_list:
            self.logger.info(f"Sample ID: {sample_id}")
            base_logger.print_space(self.logger)
            self._get_files_path(sample_id)
            self._check_files_path()
            self.sample_data[sample_id] = OneSampleData(*self.files_path)

    def _add_unspecified_sample_id_list(self):
        uniq_dir, suffix =self.import_info[0]
        self.logger.info(f"Searching sample IDs: prefix with the suffix '{suffix}' in the directory: {uniq_dir}.")
        file_list = os.listdir(uniq_dir)
        sample_id_list = [file.replace(suffix, '') for file in file_list if file.endswith(suffix)]
        self.logger.info(f"Found {len(sample_id_list)} samples.")
        base_logger.print_space(self.logger)
        self.sample_id_list.extend(sample_id_list)

    def _add_specified_sample_id_list(self, sample_id_list: list[str]):
        for sample_id in sample_id_list:
            if sample_id not in self.sample_id_list:
                self.sample_id_list.append(sample_id)
            else:
                self.logger.warning(f"WARNING: Sample ID '{sample_id}' already exists in the current instance. Skipping import.")
        base_logger.print_space(self.logger)

    def _get_files_path(self, sample_id: str):
        self.files_path = [os.path.join(in_dir, f"{sample_id}{suffix}") for in_dir, suffix in self.import_info]

    def _check_files_path(self) -> None:
        for file_path in self.files_path:
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"File does not exist: {file_path}.")

    @base_logger.prog_log(prog_name="Import sample metadata")
    def _read_sample_metadata(self, sample_metadata_path: str) -> None:
        self.sample_metadata = {}
        df = pd.read_csv(sample_metadata_path, index_col="Sample_id")
        for sample_id in self.sample_id_list:
            if sample_id not in df.index:
                self.logger.warning(f"WARNING: Sample ID '{sample_id}' not found in the sample metadata table.")
            else:
                self.sample_metadata[sample_id] = df.loc[df.index == sample_id].to_dict("records")[0]

    @base_logger.prog_log(prog_name="Pickle sample data")
    def pickle_data(self,
            save_dir: str,
            save_prefix: str = f"eDNA_samples_{date.today()}",
            overwrite: bool = False
        ):
        """
        Pickle the current sample data to a specified directory.

        :param save_dir: The directory to save the SamplesContainer instance.
        :param save_prefix: The prefix for the save .pkl file. Defaults is 'eDNA_samples_<current_date>'.
        :param overwrite: If True, overwrite the existing file. Defaults to False.
        """
        os.makedirs(save_dir, exist_ok=True)
        pickle_path = os.path.join(save_dir, f"{save_prefix}.pkl")
        if os.path.exists(pickle_path) and not overwrite:
            self.logger.warning(f"WARNING: File already exists: {pickle_path}. Data didn't saved.")
            return
        self._pickle_instance(pickle_path)

    def _pickle_instance(self, pickle_path) -> None:
        with open(pickle_path, 'wb') as f:
            pickle.dump(self, f)

    @base_logger.prog_log(prog_name="Unpickle sample data")
    def unpickle_data(self, pickle_path) -> None:
        """
        Unpickle sample data from a specified path.
    
        :param pickle_path: The path to unpickle a pre-built SamplesData pickle file.
        """
        with open(pickle_path,'rb') as file:
            self.__dict__ = pickle.load(file).__dict__

    @base_logger.prog_log(prog_name="Merge input object(s) into the current object")
    def merge_data(self, *sample_data: object) -> None:
        """
        Merge sample data from another SamplesContainer instance into the current instance.

        :param sample_data: Variable number of SampleData instances to merge from.
        """
        for data_object in sample_data:
            self._validate_input_object(data_object)
            self._merge_single_data_object(data_object)
    
    def _validate_input_object(self, data_object: object) -> None:
        if not isinstance(data_object, SampleData):
            raise TypeError("FAIL: Input object must be an instance of SamplesData.")
    
    def _merge_single_data_object(self, data_object: 'SampleData') -> None:
        new_sample_ids = [sample_id for sample_id in data_object.sample_id_list if sample_id not in self.sample_id_list]
    
        for sample_id in data_object.sample_id_list:
            if sample_id not in new_sample_ids:
                self.logger.warning(f"WARNING: Sample ID '{sample_id}' already exists in the current instance. Skipping merge.")
                continue
            self._merge_sample(sample_id, data_object)
            self._merge_metadata(sample_id, data_object)
            self.sample_id_list.append(sample_id)
    
    def _merge_sample(self, sample_id: str, data_object: 'SampleData') -> None:
        self.sample_data[sample_id] = data_object.sample_data[sample_id]
    
    def _merge_metadata(self, sample_id: str, data_object: 'SampleData') -> None:
        if not hasattr(data_object, "sample_metadata"):
            return
        if not hasattr(self, "sample_metadata"):
            self.sample_metadata = {}
        self.sample_metadata[sample_id] = data_object.sample_metadata[sample_id]