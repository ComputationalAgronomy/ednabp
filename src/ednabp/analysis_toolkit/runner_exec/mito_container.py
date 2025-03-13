import os
import warnings
from typing import override

import pandas as pd

from ..runner_build import base_logger
from .data_container import SampleData


class OneMitoData:
    """A class to handle MitoFish outputs from Excel sheets.

    Attributes:
        xls (pd.ExcelFile): Excel file object containing mitochondrial data.
        spc_smpdata_df (pd.DataFrame): DataFrame containing sample details from 'List of Sample Details' sheet.
        spc_metadata_df (pd.DataFrame): DataFrame containing species comparison data from 'Comparison of Samples' sheet.
        hap_seq (dict): Dictionary mapping haploid IDs to their corresponding sequences.
        hap_size (dict): Dictionary mapping haploid IDs to their sequence sizes.
        hap2level (dict): Dictionary mapping haploid IDs to their taxonomic information (species, genus, class, order, family).
        spc_info (dict): Dictionary containing species-specific information including water area, habitat, depth ranges,
                        IUCN status, and fisheries importance.
    """

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

    def __init__(self, mito_xlsx):
        try:
            self.xls = pd.ExcelFile(mito_xlsx)
            self._validate_sheets()
        except FileNotFoundError:
            raise FileNotFoundError(f"Excel file not found: {mito_xlsx}")

        self.read_sheets()
        self.load_data()

    def _validate_sheets(self):
        required_sheets = {
            self.SAMPLE_DETAILS_SHEET,
            self.SPECIES_COMPARISON_SHEET,
        }
        missing_sheets = required_sheets - set(self.xls.sheet_names)
        if missing_sheets:
            raise ValueError(f"Missing required sheets: {missing_sheets}")

    def read_sheets(self):
        try:
            self.spc_smpdata_df = pd.read_excel(
                self.xls, self.SAMPLE_DETAILS_SHEET
            )
            self._validate_sample_data_columns()
            self.spc_smpdata_df["Species"] = self.spc_smpdata_df[
                "Species"
            ].ffill()

            self.spc_metadata_df = pd.read_excel(
                self.xls, self.SPECIES_COMPARISON_SHEET
            )
            self._validate_metadata_columns()
            self.spc_metadata_df = self.spc_metadata_df.dropna(how="all")
            self.spc_metadata_df = self.spc_metadata_df[
                self.spc_metadata_df["Class"]
                != "Followings are non-fish species"
            ]
            self.spc_metadata_df = self.spc_metadata_df.fillna("Not record")
        except Exception as e:
            raise ValueError(f"Error reading Excel sheets: {str(e)}")

    def _validate_sample_data_columns(self):
        required_columns = {"Haploid ID", "Sequence", "Size", "Species"}
        self._check_missing_columns(self.spc_smpdata_df, required_columns)

    def _validate_metadata_columns(self):
        required_columns = set(
            self.TAXONOMIC_COLUMNS + self.SPECIES_INFO_COLUMNS
        )
        self._check_missing_columns(self.spc_metadata_df, required_columns)

    @staticmethod
    def _check_missing_columns(df: pd.DataFrame, required_columns: set):
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

    def load_data(self):
        if self.spc_smpdata_df is None or self.spc_metadata_df is None:
            raise ValueError("Must call read_sheets() before load_data()")
        self._create_haploid_mappings()
        self._create_taxonomic_mappings()
        self._create_species_info_mapping()

    def _create_haploid_mappings(self):
        hap_id = self.spc_smpdata_df["Haploid ID"]
        self.hap_seq = dict(
            zip(hap_id, self.spc_smpdata_df["Sequence"], strict=False)
        )
        self.hap_size = dict(
            zip(hap_id, self.spc_smpdata_df["Size"], strict=False)
        )

    def _create_taxonomic_mappings(self):
        hap_id = self.spc_smpdata_df["Haploid ID"]
        species = self.spc_smpdata_df["Species"]

        # Create initial hap2level with species and genus
        self.hap2level = {
            hap: {"species": spc, "genus": spc.split()[0]}
            for hap, spc in zip(hap_id, species, strict=False)
        }

        # Add additional taxonomic information
        spc2level = (
            self.spc_metadata_df.filter(items=self.TAXONOMIC_COLUMNS)
            .set_index("Scientific Name")
            .T.to_dict("dict")
        )

        for hap, level in self.hap2level.items():
            self.hap2level[hap].update(spc2level[level["species"]])
            self.hap2level[hap] = {
                k.lower(): v for k, v in self.hap2level[hap].items()
            }

    def _create_species_info_mapping(self):
        self.spc_info = (
            self.spc_metadata_df.filter(items=self.SPECIES_INFO_COLUMNS)
            .set_index("Scientific Name")
            .T.to_dict("dict")
        )


class MitoData(SampleData):
    """
    A class extends SampleData to provide specific functionality for
    processing MitoFish outputs stored in Excel format.
    """

    def __init__(self, no_verbose=False):
        super().__init__(no_verbose)

    @override
    def import_data(
        self,
        mito_xlsx_dir: str,
        mito_xlsx_suffix: str = ".xlsx",
        sample_id_list: list[str] | None = None,
        sample_metadata_path: str | None = None,
    ):
        self.import_sample_id_list = []
        self.in_dir = mito_xlsx_dir
        self.in_suffix = mito_xlsx_suffix

        self._read_sample_data(sample_id_list)

        if sample_metadata_path is not None:
            self._read_sample_metadata(sample_metadata_path)

        self.sample_id_list.extend(self.import_sample_id_list)

    @override
    @base_logger.prog_log(prog_name="Import data from MitoFish outputs")
    def _read_sample_data(self, sample_id_list: list[str] | None):
        base_logger.print_space(self.logger)
        if sample_id_list is None:
            self.logger.info("No sample id list provided.")
            self._add_unspecified_sample_id_list()
        else:
            self.logger.info("Specified sample id list.")
            self._add_specified_sample_id_list(sample_id_list)

        for sample_id in self.import_sample_id_list:
            self.logger.info(f"Sample ID: {sample_id}")
            base_logger.print_space(self.logger)
            one_mito_data = OneMitoData(self._get_file_path(sample_id))
            self.spc_info.update(one_mito_data.spc_info)
            delattr(one_mito_data, "spc_info")
            self.sample_data[sample_id] = one_mito_data

    @override
    def _add_unspecified_sample_id_list(self):
        self.logger.info(
            f"Searching sample IDs: prefix with the suffix '{self.in_suffix}' in the directory: {self.in_dir}."
        )
        sample_id_list = [
            os.path.splitext(file)[0]
            for file in os.listdir(self.in_dir)
            if file.endswith(self.in_suffix)
        ]
        self.logger.info(f"Found {len(sample_id_list)} samples.")
        for sample_id in sample_id_list:
            if sample_id in self.sample_id_list:
                self.logger.warning(
                    f"WARNING: Sample ID '{sample_id}' already exists in the current instance. "
                    "Skipping import."
                )
                continue
            self.import_sample_id_list.append(sample_id)
        base_logger.print_space(self.logger)

    @override
    def _add_specified_sample_id_list(self, sample_id_list: list[str]):
        for sample_id in sample_id_list:
            if not self._is_valid_sample_file(sample_id):
                self.logger.warning(
                    f"WARNING: Sample ID '{sample_id}' (file: '{sample_id + self.in_suffix}') "
                    f"not found in the directory. Skipping import."
                )
                continue

            if sample_id in self.sample_id_list:
                self.logger.warning(
                    f"WARNING: Sample ID '{sample_id}' already exists in the current instance. "
                    "Skipping import."
                )
                continue

            self.import_sample_id_list.append(sample_id)
        base_logger.print_space(self.logger)

    def _get_file_path(self, sample_id: str) -> str:
        return os.path.join(self.in_dir, f"{sample_id}{self.in_suffix}")

    def _is_valid_sample_file(self, sample_id: str):
        return os.path.exists(self._get_file_path(sample_id))
