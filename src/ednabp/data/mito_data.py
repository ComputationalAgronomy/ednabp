import os

import pandas as pd

from ..common import base_logger
from .bp_data import BPData

SAMPLE_DETAILS_SHEET = "List of Sample Details"
SPECIES_COMPARISON_SHEET = "Comparison of Samples"

TAXONOMIC_COLUMNS = ["Class", "Order", "Family", "Scientific Name"]
SPECIES_INFO_COLUMNS = [
    "Scientific Name",
    "Water area",
    "Habitat",
    "DepthS",
    "DepthD",
    "IUCN Red List Status",
    "Importance in Fisheries",
]


class SingleMitoData:
    """
    Attributes:
        hap_seq (dict): Dictionary mapping haploid IDs to their corresponding sequences.
        hap_size (dict): Dictionary mapping haploid IDs to their sequence sizes.
        hap2level (dict): Dictionary mapping haploid IDs to their taxonomic information (species, genus, class, order, family).
    """

    pass


class MitoFileReader:
    """A class to handle MitoFish outputs from Excel sheets.

    Attributes:
        xls (pd.ExcelFile): Excel file object containing mitochondrial data.
        spc_smpdata_df (pd.DataFrame): DataFrame containing sample details from 'List of Sample Details' sheet.
        spc_smpdata_df_dict (Dictionary): Dictionary containing sample name specific details from spc_smpdata_df.
        spc_metadata_df (pd.DataFrame): DataFrame containing species comparison data from 'Comparison of Samples' sheet.
        spc_info (dict): Dictionary containing species-specific information including water area, habitat, depth ranges,
                        IUCN status, and fisheries importance.
    """

    def __init__(self, mito_xlsx):
        try:
            self.xls = pd.ExcelFile(mito_xlsx)
            self._validate_sheets()
        except FileNotFoundError:
            raise FileNotFoundError(f"Excel file not found: {mito_xlsx}")
        self.mito_data = {}

        self.read_sheets()
        self.load_data()

    def _validate_sheets(self):
        required_sheets = {
            SAMPLE_DETAILS_SHEET,
            SPECIES_COMPARISON_SHEET,
        }
        missing_sheets = required_sheets - set(self.xls.sheet_names)
        if missing_sheets:
            raise ValueError(f"Missing required sheets: {missing_sheets}")

    def read_sheets(self):
        try:
            self.spc_smpdata_df = pd.read_excel(self.xls, SAMPLE_DETAILS_SHEET)
            self.validate_sample_data_columns()
            self.spc_smpdata_df = self.spc_smpdata_df[
                ~self.spc_smpdata_df["Species"].isin(["Species", "Skip"])
            ]
            self.spc_smpdata_df["Sample name"] = self.spc_smpdata_df[
                "Sample name"
            ].ffill()
            self.spc_smpdata_df["Species"] = self.spc_smpdata_df[
                "Species"
            ].ffill()
            self.spc_smpdata_df_dict = {
                sample_name: self.spc_smpdata_df[
                    self.spc_smpdata_df["Sample name"] == sample_name
                ]
                for sample_name in self.spc_smpdata_df["Sample name"].unique()
            }

            self.spc_metadata_df = pd.read_excel(
                self.xls, SPECIES_COMPARISON_SHEET
            )
            self.validate_metadata_columns()
            self.spc_metadata_df = self.spc_metadata_df.dropna(how="all")
            self.spc_metadata_df = self.spc_metadata_df[
                self.spc_metadata_df["Class"]
                != "Followings are non-fish species"
            ]
            self.spc_metadata_df = self.spc_metadata_df.fillna("Not record")
        except Exception as e:
            raise ValueError(f"Error reading Excel sheets: {str(e)}")

    def validate_sample_data_columns(self):
        required_columns = {"Haploid ID", "Sequence", "Size", "Species"}
        self.check_missing_columns(self.spc_smpdata_df, required_columns)

    def validate_metadata_columns(self):
        required_columns = set(TAXONOMIC_COLUMNS + SPECIES_INFO_COLUMNS)
        self.check_missing_columns(self.spc_metadata_df, required_columns)

    @staticmethod
    def check_missing_columns(df: pd.DataFrame, required_columns: set):
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

    def load_data(self):
        if self.spc_smpdata_df is None or self.spc_metadata_df is None:
            raise ValueError("Must call read_sheets() before load_data()")
        for sample_name, df in self.spc_smpdata_df_dict.items():
            one_sample_data = SingleMitoData()
            one_sample_data.hap_seq = self.haploid2sequnce_mappings(df)
            one_sample_data.hap_size = self.haploid2size_mappings(df)
            one_sample_data.hap2level = self.taxonomic_mappings(df)
            self.mito_data[sample_name] = one_sample_data

        self.create_species_info_mapping()

    def haploid2sequnce_mappings(self, df):
        hap_seq = dict(zip(df["Haploid ID"], df["Sequence"], strict=False))
        return hap_seq

    def haploid2size_mappings(self, df):
        hap_size = dict(zip(df["Haploid ID"], df["Size"], strict=False))
        return hap_size

    def taxonomic_mappings(self, df):
        hap_id = df["Haploid ID"]
        species = df["Species"]

        # Create initial hap2level with species and genus
        hap2level = {
            hap: {"species": spc, "genus": spc.split()[0]}
            for hap, spc in zip(hap_id, species, strict=False)
        }

        # Add additional taxonomic information
        spc2level = (
            self.spc_metadata_df.filter(items=TAXONOMIC_COLUMNS)
            .set_index("Scientific Name")
            .T.to_dict("dict")
        )

        for hap, level in hap2level.items():
            hap2level[hap].update(spc2level[level["species"]])
            hap2level[hap] = {k.lower(): v for k, v in hap2level[hap].items()}

        return hap2level

    def create_species_info_mapping(self):
        self.spc_info = (
            self.spc_metadata_df.filter(items=SPECIES_INFO_COLUMNS)
            .set_index("Scientific Name")
            .T.to_dict("dict")
        )


class MitoData(BPData):
    """
    A class extends SampleData to provide specific functionality for
    processing MitoFish outputs stored in Excel format.
    """

    def __init__(self, verbose: bool = True):
        super().__init__(verbose)

    def import_data(
        self,
        mito_xlsx_dir: str,
        mito_xlsx_suffix: str = ".xlsx",
        file_list: list[str] | None = None,
        sample_metadata_path: str | None = None,
        date_column: str = "date",
        date_format: str = "%Y-%m",
    ):
        self.import_file_list = []
        self.import_sample_id_list = []
        self.in_dir = mito_xlsx_dir
        self.in_suffix = mito_xlsx_suffix

        self.read_sample_data(file_list)

        if sample_metadata_path is not None:
            self.import_metadata(
                sample_metadata_path, date_column, date_format
            )

        self.sample_id_list.extend(self.import_sample_id_list)

    @base_logger.prog_log(prog_name="Import data from MitoFish outputs")
    def read_sample_data(self, file_list: list[str] | None):
        if file_list is None:
            self.logger.info("No sample id list provided.")
            self.add_unspecified_file_list()
        else:
            self.logger.info("Specified sample id list.")
            self.add_specified_file_list(file_list)

        for file_name in self.import_file_list:
            self.logger.info(f"Sample ID: {file_name}")
            read_results = MitoFileReader(self.get_file_path(file_name))
            self.logger.info(f"Read {len(read_results.mito_data)} sample(s)")
            for sample_id in read_results.mito_data.keys().copy():
                if (
                    sample_id in self.sample_id_list
                    or sample_id in self.import_sample_id_list
                ):
                    self.logger.warning(
                        f"WARNING: Sample ID '{sample_id}' in the file {file_name} already exists in the current instance. "
                        "Skipping import."
                    )
                    read_results.mito_data.pop(sample_id, None)
                else:
                    self.import_sample_id_list.append(sample_id)

            self.spc_info.update(read_results.spc_info)
            self.sample_data.update(read_results.mito_data)

    def add_unspecified_file_list(self):
        self.logger.info(
            f"Searching sample IDs: prefix with the suffix '{self.in_suffix}' in the directory: {self.in_dir}."
        )
        self.import_file_list = [
            os.path.splitext(file)[0]
            for file in os.listdir(self.in_dir)
            if file.endswith(self.in_suffix)
        ]
        self.logger.info(f"Found {len(self.import_file_list)} file(s).")

    def add_specified_file_list(self, file_list: list[str]):
        for file_name in file_list:
            if not self.is_valid_file(self.get_file_path(file_name)):
                self.logger.warning(
                    f"WARNING: File name '{file_name}' (file: '{file_name + self.in_suffix}') "
                    f"not found in the directory. Skipping import."
                )
                continue
            self.import_file_list.append(file_name)

    def get_file_path(self, sample_id: str) -> str:
        return os.path.join(self.in_dir, f"{sample_id}{self.in_suffix}")
