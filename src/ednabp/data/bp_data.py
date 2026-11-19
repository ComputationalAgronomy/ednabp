import os
import pickle
import sys
from datetime import date

import chardet
import duckdb
import pandas as pd

from ..common import base_logger

FILL_NA = "N/A"


class SingleBPData:
    """
    Container for handling ZOTU data from taxa CSV files.
    This class reads ZOTU sequences, sizes, and taxonomy information from a single taxa.csv file.

    :param taxa_csv: Path to the taxa CSV file containing ZOTU information.

    :attribute hap_seq: A dictionary containing ZOTU sequences.
    :attribute hap_size: A dictionary containing ZOTU sizes.
    :attribute hap2level: A dictionary mapping ZOTUs to taxonomic levels.
    """

    DESIRED_LEVELS = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
        "domain",
    ]

    def __init__(self, taxa_csv: str | None = None, sample_id=None):
        """
        :param taxa_csv: addhap output of one sample.
        :param sample_id: The sample's abundance column, used if the table
            has no "size" column (addhap names the column after the sample).
        """
        self.hap_seq = {}
        self.hap_size = {}
        self.hap2level = {}
        if taxa_csv is not None:
            self._read_csv(taxa_csv, sample_id)

    @classmethod
    def from_frame(cls, df: pd.DataFrame, size_column: str = "size"):
        """
        Build from a taxa table already in memory; `size_column` holds the
        abundances (a sample's column of a pooled table).
        """
        data = cls()
        data._load_frame(df, size_column)
        return data

    @base_logger.prog_log(prog_name="Read Taxa CSV")
    def _read_csv(self, taxa_csv: str, sample_id=None):
        try:
            df = pd.read_csv(taxa_csv)
            if df.empty:
                base_logger.logger.warning("Taxa CSV is empty")
                return
            df.columns = [str(c) for c in df.columns]
            size_column = "size"
            if "size" not in df.columns and str(sample_id) in df.columns:
                size_column = str(sample_id)
            self._load_frame(df, size_column)
        except (pd.errors.EmptyDataError, FileNotFoundError) as e:
            base_logger.logger.warning(f"Error reading taxa CSV: {e}")

    def _load_frame(self, df: pd.DataFrame, size_column: str = "size"):
        # Read sequences if available
        if "sequence" in df.columns:
            self.hap_seq = dict(
                zip(df["qseqid"], df["sequence"], strict=False)
            )

        # Read sizes if available
        if size_column in df.columns:
            self.hap_size = dict(
                zip(df["qseqid"], df[size_column], strict=False)
            )

        # Read taxonomy levels
        available_levels = [
            col for col in self.DESIRED_LEVELS if col in df.columns
        ]
        if available_levels:
            df_indexed = df.set_index("qseqid")[available_levels]
            self.hap2level = df_indexed.to_dict("index")


class BPData:
    """
    A class for managing sample data storage from taxa CSV files.
    It provides methods for importing, pickling, unpickling, merging.

    :attribute sample_data: A dictionary to store sample data, the key is sample ID, and the value is a SingleBPData instance.
    :attribute sample_metadata: A dictionary to store sample information, the key is sample ID, and the value is a dictionary containing sample metadata.
    :attribute sample_id_list: A list to store sample IDs.
    :attribute verbose: A boolean flag to control logging verbosity. Default is True.
    """

    SAMPLE_ID_COLUMN = "sample_id"

    def __init__(self, verbose=True):
        self.sample_data = {}
        self.sample_metadata = {}
        self.spc_info = {}
        self.import_info = []
        self.sample_id_list = []
        self.verbose = verbose
        self.logger = base_logger.logger

        if self.verbose:
            self.logger.setLevel("INFO")
        else:
            self.logger.setLevel("WARNING")

    def import_data(
        self,
        results_dir: str,
        sample_id_list: list[str] | None = None,
        file_suffix: str = "_blast.csv",
    ):
        """
            Import sample data from taxa CSV files.

            :param results_dir: Parent directory containing output CSV files.
            :param sample_id_list: List of sample IDs to import. If None, imports all available.
        :param file_suffix: File suffix for blast files. Default is '_blast.csv'.
        """
        self.import_sample_id_list = []
        self.results_dir = results_dir

        self.file_suffix = file_suffix
        self.read_sample_data(sample_id_list)
        self.sample_id_list.extend(self.import_sample_id_list)

    @base_logger.prog_log(prog_name="Import pooled OTU table")
    def import_otu_table(
        self,
        otu_table_path: str,
        sample_id_list: list[str] | None = None,
        min_count: int = 1,
    ):
        """
        Import samples from one pooled taxa table (the pipeline's
        `blast/_pooled.csv` when run with a pool stage). Each row is an
        OTU; the columns after "sequence" are the samples' abundances.

        Each sample keeps only the OTUs with at least `min_count` reads in
        it, so that occurrence and richness count only OTUs present in the
        sample. Raise `min_count` to drop low counts from cross-sample
        leakage (tag jumping).

        :param otu_table_path: Path to the pooled taxa table.
        :param sample_id_list: Samples (columns) to import. If None, imports all.
        :param min_count: Minimum abundance of an OTU in a sample. Default: 1.
        """
        df = pd.read_csv(otu_table_path)
        if "sequence" not in df.columns:
            raise ValueError(
                f"{otu_table_path} has no 'sequence' column; the sample"
                " columns follow it."
            )
        sample_columns = [
            str(c)
            for c in df.columns[list(df.columns).index("sequence") + 1 :]
        ]
        df.columns = [str(c) for c in df.columns]

        requested = (
            sample_columns if sample_id_list is None else sample_id_list
        )
        self.import_sample_id_list = []
        self.results_dir = os.path.dirname(otu_table_path)
        for sample_id in requested:
            sample_id = str(sample_id)
            if sample_id in self.sample_id_list:
                self.logger.warning(
                    f"Sample ID '{sample_id}' already exists in the current instance. "
                    "Skipping import."
                )
                continue
            if sample_id not in sample_columns:
                self.logger.warning(
                    f"Sample ID '{sample_id}' is not a column of {otu_table_path}. "
                    "Skipping import."
                )
                continue
            present = df[df[sample_id] >= max(min_count, 1)]
            self.sample_data[sample_id] = SingleBPData.from_frame(
                present, size_column=sample_id
            )
            self.import_sample_id_list.append(sample_id)

        self.logger.info(
            f"Imported {len(self.import_sample_id_list)} samples from"
            f" {otu_table_path}."
        )
        self.sample_id_list.extend(self.import_sample_id_list)

    @base_logger.prog_log(prog_name="Import sample sequence data")
    def read_sample_data(self, sample_id_list: list[str] | None):
        if sample_id_list is None:
            self.logger.info("No sample id list provided.")
            self.add_unspecified_sample_id_list()
        else:
            self.logger.info("Specified sample id list.")
            self.add_specified_sample_id_list(sample_id_list)

        for sample_id in self.import_sample_id_list:
            self.logger.info(f"Sample ID: {sample_id}")
            taxa_file = os.path.join(
                self.results_dir, f"{sample_id}{self.file_suffix}"
            )
            self.sample_data[sample_id] = SingleBPData(taxa_file, sample_id)

    def add_unspecified_sample_id_list(self):
        self.logger.info(
            f"Searching sample IDs with suffix '{self.file_suffix}' in directory: {self.results_dir}."
        )
        file_list = os.listdir(self.results_dir)
        sample_id_list = [
            file.replace(self.file_suffix, "")
            for file in file_list
            if file.endswith(self.file_suffix)
        ]
        self.logger.info(f"Found {len(sample_id_list)} samples.")
        for sample_id in sample_id_list:
            if sample_id in self.sample_id_list:
                self.logger.warning(
                    f"Sample ID '{sample_id}' already exists in the current instance. "
                    "Skipping import."
                )
                continue
            taxa_file = os.path.join(
                self.results_dir, f"{sample_id}{self.file_suffix}"
            )
            if self.is_valid_file(taxa_file):
                self.import_sample_id_list.append(sample_id)
            else:
                self.logger.warning(
                    f"File {taxa_file} doesn't exist. Skipping import."
                )

    def add_specified_sample_id_list(self, sample_id_list: list[str]):
        for sample_id in sample_id_list:
            if sample_id in self.sample_id_list:
                self.logger.warning(
                    f"Sample ID '{sample_id}' already exists in the current instance. "
                    "Skipping import."
                )
                continue
            taxa_file = os.path.join(
                self.results_dir, f"{sample_id}{self.file_suffix}"
            )
            if self.is_valid_file(taxa_file):
                self.import_sample_id_list.append(sample_id)
            else:
                self.logger.warning(
                    f"File {taxa_file} doesn't exist. Skipping import."
                )

    def is_valid_file(self, file_path: str):
        return os.path.exists(file_path)

    @base_logger.prog_log(prog_name="Import sample metadata")
    def import_metadata(
        self,
        sample_metadata_path: str | None = None,
    ) -> None:
        with open(sample_metadata_path, "rb") as f:
            result = chardet.detect(f.read())
        df = pd.read_csv(
            sample_metadata_path,
            index_col=self.SAMPLE_ID_COLUMN,
            encoding=result["encoding"],
        )

        for sample_id in self.import_sample_id_list:
            if str(sample_id) not in df.index.astype(str):
                self.logger.warning(
                    f"Sample ID '{sample_id}' not found in the sample metadata table."
                )
            else:
                self.sample_metadata[sample_id] = df.loc[
                    df.index.astype(str) == str(sample_id)
                ].to_dict("records")[0]

    @base_logger.prog_log(prog_name="Import species information")
    def import_spc_info(self, fishbase_db_path, stock_db_path):
        spc_info = {}
        conn = duckdb.connect()
        try:
            link_fishbase = conn.from_parquet(fishbase_db_path)
            link_stock = conn.from_parquet(stock_db_path)
            IUCN_levels = {
                "N.E.": "Not Evaluated",
                "DD": "Data Deficient",
                "N.A.": "Not Available",
                "LC": "Least Concern",
                "LR/lc": "Lower Risk: least concern",
                "NT": "Near Threatened",
                "LR/cd": "Vulnerable",
                "VU": "Vulnerable",
                "LR/nt": "Lower Risk: near threatened",
                "EN": "Endangered",
                "EX": "Extinct",
                "EW": "Extinct in the Wild",
                "CR": "Critically Endangered",
            }

            all_species = {
                level["species"]
                for sample_id in self.sample_id_list
                for level in self.sample_data[sample_id].hap2level.values()
            }

            for species_name in all_species:
                if not isinstance(species_name, str):
                    continue
                (genus_name, *species_subnames) = species_name.split("_")
                species_subname = (
                    species_subnames[0] if len(species_subnames) > 0 else "sp"
                )
                species_subname = species_subname.replace("'", "").replace(
                    "-", "_"
                )
                fb_data = (
                    link_fishbase.filter(
                        f"Genus = '{genus_name}' AND Species = '{species_subname}'"
                    )
                    .project(
                        "Fresh, Brack, Saltwater, DemersPelag, \
                    DepthRangeShallow, DepthRangeDeep, Importance, SpecCode"
                    )
                    .fetchone()
                )

                if fb_data is not None:
                    water = ""
                    if fb_data[0] == 1:
                        water += "Fresh Water; "
                    if fb_data[1] == 1:
                        water += "Salt Water; "
                    if fb_data[2] == 1:
                        water += "Brack Water; "
                    habitat = fb_data[3] if fb_data[3] is not None else FILL_NA
                    depth_s = fb_data[4] if fb_data[4] is not None else FILL_NA
                    depth_d = fb_data[5] if fb_data[5] is not None else FILL_NA
                    fisheries_role = (
                        fb_data[6] if fb_data[6] is not None else FILL_NA
                    )
                    iucn_data = (
                        link_stock.filter(f"SpecCode = {fb_data[7]}")
                        .project("IUCN_Code")
                        .fetchone()
                    )
                    iucn_lv = (
                        IUCN_levels[iucn_data[0]]
                        if iucn_data[0] in IUCN_levels
                        else FILL_NA
                    )
                else:
                    (
                        water,
                        habitat,
                        depth_s,
                        depth_d,
                        fisheries_role,
                        iucn_lv,
                    ) = [FILL_NA] * 6
                spc_info[species_name] = {
                    "Water area": water,
                    "Habitat": habitat,
                    "DepthS": depth_s,
                    "DepthD": depth_d,
                    "Importance in Fisheries": fisheries_role,
                    "IUCN Red List Status": iucn_lv,
                }

        finally:
            conn.close()

        self.spc_info.update(spc_info)

    @base_logger.prog_log(prog_name="Pickle sample data")
    def pickle_data(
        self,
        save_dir: str,
        save_prefix: str = f"eDNA_samples_{date.today()}",
        overwrite: bool = False,
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
            self.logger.warning(
                f"File already exists: {pickle_path}. Data didn't saved."
            )
            return
        self.pickle_instance(pickle_path)

    def pickle_instance(self, pickle_path) -> None:
        with open(pickle_path, "wb") as f:
            pickle.dump(self, f)

    @base_logger.prog_log(prog_name="Unpickle sample data")
    def unpickle_data(self, pickle_path) -> None:
        """
        Unpickle sample data from a specified path.

        :param pickle_path: The path to unpickle a pre-built SamplesData pickle file.
        """
        with open(pickle_path, "rb") as file:
            self.__dict__ = pickle.load(file).__dict__

    @base_logger.prog_log(
        prog_name="Merge input object(s) into the current object"
    )
    def merge_data(self, *sample_data: object) -> None:
        """
        Merge sample data from another SamplesContainer instance into the current instance.

        :param sample_data: Variable number of SampleData instances to merge from.
        """
        for data_object in sample_data:
            self.validate_input_object(data_object)
            self.merge_single_data_object(data_object)

    def validate_input_object(self, data_object: object) -> None:
        if not isinstance(data_object, BPData):
            self.logger.error(
                "FAIL: Input object must be an instance of BPData."
            )
            sys.exit(1)

    def merge_single_data_object(self, data_object: "BPData") -> None:
        new_sample_ids = [
            sample_id
            for sample_id in data_object.sample_id_list
            if sample_id not in self.sample_id_list
        ]

        for sample_id in data_object.sample_id_list:
            if sample_id not in new_sample_ids:
                self.logger.warning(
                    f"Sample ID '{sample_id}' already exists in the current instance. Skipping merge."
                )
                continue
            self.merge_sample(sample_id, data_object)
            self.merge_metadata(sample_id, data_object)
            self.sample_id_list.append(sample_id)
        self.merge_spc_info(data_object)

    def merge_sample(self, sample_id: str, data_object: "BPData") -> None:
        self.sample_data[sample_id] = data_object.sample_data[sample_id]

    def merge_metadata(self, sample_id: str, data_object: "BPData") -> None:
        if not hasattr(data_object, "sample_metadata"):
            return
        self.sample_metadata[sample_id] = data_object.sample_metadata[
            sample_id
        ]

    def merge_spc_info(self, data_object: "BPData") -> None:
        if not hasattr(data_object, "spc_info"):
            return
        self.spc_info.update(data_object.spc_info)
