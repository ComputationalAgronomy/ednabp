import os
from collections import defaultdict
from itertools import product
from typing import Literal

import numpy as np
import pandas as pd

from ...common import base_logger, base_writer
from ..plot.base_plotter import VALUE_COLUMN

UNIT_COLUMN = "unit"
SAMPLE_ID_COLUMN = "sample_id"
FILL_NA = "N/A"


def get_unit_name(level_dict: dict, unit_level: str, hap: str) -> str:
    if unit_level in level_dict:
        return level_dict[unit_level]
    if unit_level == "haplotype":
        return f"{level_dict['species']}_{hap}"
    raise ValueError(f"Invalid unit_level: {unit_level}")


@base_logger.prog_log("Export dataframe to CSV file")
def export_df(metric_df, save_dir, file_name, overwrite, logger):
    output_path = os.path.join(save_dir, file_name)
    if os.path.exists(output_path) and not overwrite:
        logger.warning(f"File already exists: {output_path}. Stop saving.")
        return
    os.makedirs(save_dir, exist_ok=True)
    metric_df.to_csv(output_path, index=False)
    logger.info(f"Dataframe exported to: {output_path}")


class Writer(base_writer.BaseWriter):
    """
    An abstract class for running diversity metrics related analysis
    Metrics include:
        - abundance: sum of sequence read for each target taxon. (e.g. How many sequences of fish got in sampleA?)
        - richness: number of unique units for each target taxon. (e.g. How many different fish species got in sampleA?)
    """

    def __init__(self, data, verbose=False):
        super().__init__(data, verbose)

    @base_logger.prog_log("Write richness table")
    def richness(
        self,
        taxa_lv: str,
        unit_lv: str = "species",
        save_dir: str = None,
        overwrite: bool = False,
        sample_id_list: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Write the richness data to a CSV file.

        :param save_dir: The directory path to save the CSV file.
        :param taxa_level: The taxa level to group units (e.g., species, family, etc.).
        :param unit_level: The taxa level to calculate richness (e.g., species, haplotype, etc.). Default is "species".
        :param sample_id_list: List of sample IDs to import. Default: None (all available sample IDs will be used).
        """
        self.taxa_lv = taxa_lv
        self.unit_lv = unit_lv
        self.occur_df = pd.DataFrame()
        self.richness_df = pd.DataFrame()
        self.load_sample_id_list(sample_id_list)
        self.create_richness_df()
        if save_dir is not None:
            export_df(
                self.richness_df,
                save_dir,
                f"{taxa_lv}_{unit_lv}_richness.csv",
                overwrite,
                self.config.logger,
            )
        return self.richness_df

    @base_logger.prog_log("Write Abundance table")
    def abundance(
        self,
        taxa_lv: str,
        process: Literal["norm", "log"] | None = None,
        save_dir: str = None,
        overwrite: bool = False,
        sample_id_list: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Write the abundance data to a CSV file.

        :param save_dir: The directory to save the CSV file.
        :param taxa_level: The name of the level to target (e.g., species, family, etc.).
        :param process: Optional data processing method:
            - 'norm': Normalize the abundance values
            - 'log': Apply log transformation
            - None: No processing (default)
        :param sample_id_list: A list of sample IDs to write.
        """
        self.taxa_lv = taxa_lv
        self.process = process
        self.abundance_df = pd.DataFrame()
        self.load_sample_id_list(sample_id_list)
        self.create_abundance_df()
        if save_dir is not None:
            process = "raw" if process is None else f"{process}"
            export_df(
                self.abundance_df,
                save_dir,
                f"{taxa_lv}_{process}_abundance.csv",
                overwrite,
                self.config.logger,
            )
        return self.abundance_df

    @base_logger.prog_log("Write detection probability table")
    def detectprob(
        self,
        taxa_lv: str,
        detectprob_columns: str | list[str] = "sample",
        save_dir: str = None,
        overwrite: bool = False,
        sample_id_list: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Write the detect probability data to a CSV file.

        :param save_dir: The directory to save the CSV file.
        :param taxa_level: The name of the level to target (e.g., species, family, etc.).
        :param sample_column: The column name(s) to calculate detection probability in the metadata csv file
        :param sample_id_list: A list of sample IDs to write.
        """
        self.taxa_lv = taxa_lv
        self.taxa_occur = pd.DataFrame()
        self.dp_df = pd.DataFrame()
        self.add_dp_col(detectprob_columns)

        self.load_sample_id_list(sample_id_list)
        self.create_dp_df()
        if save_dir is not None:
            export_df(
                self.dp_df,
                save_dir,
                f"{taxa_lv}_detectprob_{'_'.join(detectprob_columns)}.csv",
                overwrite,
                self.config.logger,
            )
        return self.dp_df

    @base_logger.prog_log("Calculate taxa richness and create dataframe")
    def create_richness_df(self):
        for sample_id in self.sample_id_used:
            self.config.logger.info(f"Sample ID: {sample_id}")
            self.get_sample_units_occur(sample_id, self.unit_lv)
            self.update_metric_df(
                sample_id,
                self.units_occur,
                [self.taxa_lv, UNIT_COLUMN, VALUE_COLUMN],
                "occur_df",
            )
        # self._filter_richness_df()
        self.convert_units_occur_to_taxa_richness()  # columns: taxa_level, Richness, Sample_id
        self.add_sample_metadata(
            "richness_df"
        )  # columns: taxa_level, Richness, Sample_id, Site, Year, Month, Sample
        self.add_spc_info("richness_df")

    @base_logger.prog_log("Calculate taxa abundance and create dataframe")
    def create_abundance_df(self):
        for sample_id in self.sample_id_used:
            self.config.logger.info(f"Sample ID: {sample_id}")
            self.get_sample_taxa_abundance(sample_id)
            match self.process:
                case None:
                    pass
                case "norm":
                    self.normalize_taxa_abundance()
                case "log":
                    self.log_taxa_abundance()
                case _:
                    raise ValueError(f"Invalid process: {self.process}")

            self.update_metric_df(
                sample_id,
                self.taxa_abundance.items(),
                [self.taxa_lv, VALUE_COLUMN],
                "abundance_df",
            )  # columns: taxa_level, Abundance, Sample_id
        self.add_sample_metadata(
            "abundance_df"
        )  # columns: taxa_level, Abundance, Sample_id, Site, Year, Month, Sample
        self.add_spc_info("abundance_df")

    @base_logger.prog_log(
        "Calculate taxa detection probability and create dataframe"
    )
    def add_dp_col(self, dp_col):
        if type(dp_col) is list:
            self.dp_col = dp_col
        elif type(dp_col) in (str, int):
            self.dp_col = [dp_col]
        else:
            raise ValueError(f"Invalid detectprob_column: {dp_col}")

    def create_dp_df(self):
        for sample_id in self.sample_id_used:
            self.config.logger.info(f"Sample ID: {sample_id}")
            self.get_sample_units_occur(sample_id, self.taxa_lv)
            self.update_metric_df(
                sample_id,
                self.units_occur,
                [self.taxa_lv, UNIT_COLUMN, VALUE_COLUMN],
                "taxa_occur",
            )  # columns: taxa_level, unit, detect_prob, Sample_id
        self.fill_non_detect_zero(self.taxa_lv)
        self.add_sample_metadata(
            "taxa_occur"
        )  # columns: taxa_level, detect_prob, Sample_id, Site, Year, Month, Sample
        self.convert_taxa_occur_to_taxa_dp()  # columns: taxa_level, detect_prob, Site, Year, Month
        self.add_spc_info("dp_df")

    def get_sample_units_occur(self, sample_id: str, unit_lv: str):
        """
        Get the appear units for a given sample.
        """
        self.units_occur = []
        for hap, lv_dict in self.data.sample_data[sample_id].hap2level.items():
            if self.taxa_lv not in lv_dict:
                raise ValueError(f"Invalid taxa_lv: {self.taxa_lv}")
            taxon_name = lv_dict[self.taxa_lv]
            unit_name = get_unit_name(lv_dict, unit_lv, hap)
            unit = (taxon_name, unit_name, 1)
            if unit not in self.units_occur:
                self.units_occur.append(
                    unit
                )  # e.g. [('FamA', 'SpcA1', 1), ('FamA', 'SpcA2', 1)]

    def update_metric_df(self, *args):
        sample_id, metric_data, metric_columns_name, metric_df_name = args
        df = pd.DataFrame(metric_data, columns=metric_columns_name)
        if sample_id is not None:
            df[SAMPLE_ID_COLUMN] = sample_id
        metric_df = getattr(self, metric_df_name)
        updated_metric_df = pd.concat([metric_df, df], ignore_index=True)
        setattr(self, metric_df_name, updated_metric_df)

    def add_sample_metadata(self, metric_df_name):
        if not hasattr(self.data, "sample_metadata"):
            self.config.logger.warning(
                "No attribute `sample_metadata` found. Skipping."
            )
            return
        metric_df = getattr(self, metric_df_name)
        metadata_df = (
            pd.DataFrame.from_dict(self.data.sample_metadata, orient="index")
            .reset_index()
            .rename(columns={"index": SAMPLE_ID_COLUMN})
        )
        updated_metric_df = pd.merge(
            left=metric_df,
            right=metadata_df,
            on=SAMPLE_ID_COLUMN,
            how="left",
        )
        if updated_metric_df.isna().any().any():
            self.config.logger.warning(
                f"Some samples are missing metadata. Filling them with '{FILL_NA}'."
            )
            updated_metric_df = updated_metric_df.fillna(FILL_NA)
        setattr(self, metric_df_name, updated_metric_df)

    def add_spc_info(self, metric_df_name):
        if self.taxa_lv != "species":
            return
        if not hasattr(self.data, "spc_info"):
            self.config.logger.warning("No species info found. Skipping.")
            return
        metric_df = getattr(self, metric_df_name)
        spc_info_df = (
            pd.DataFrame.from_dict(self.data.spc_info, orient="index")
            .reset_index()
            .rename(columns={"index": "species"})
        )
        updated_metric_df = pd.merge(
            left=metric_df, right=spc_info_df, on="species", how="left"
        )
        if updated_metric_df.isna().any().any():
            self.config.logger.warning(
                f"Some samples are missing metadata. Filling them with '{FILL_NA}'."
            )
            updated_metric_df = updated_metric_df.fillna(FILL_NA)
        setattr(self, metric_df_name, updated_metric_df)

    # def _filter_richness_df(self, site_occur_thres: int = 0, sample_occur_thres: int = 0):
    #     sp_list = self.richness_df.groupby(["Unit"])[["Site","Year","Month","Sample"]].nunique()
    #     sp_list = sp_list[sp_list["Site"]>0] # Filtering by number of sites
    #     sp_list = sp_list[sp_list["Sample"]>0] # Filtering by number of samples
    #     survive = np.array([])
    #     for m in sp_list.index:
    #         survive = np.append(survive, np.where(self.richness_df[sp_list.index.name]==m)[0])
    #     self.richness_df = self.richness_df.iloc[survive,:].reset_index(drop=True)

    def convert_units_occur_to_taxa_richness(self):
        self.richness_df = (
            self.occur_df.groupby([self.taxa_lv, SAMPLE_ID_COLUMN])[
                VALUE_COLUMN
            ]
            .sum()
            .reset_index()
        )

    def get_sample_taxa_abundance(self, sample_id: str):
        self.taxa_abundance = defaultdict(int)
        for hap, level_dict in self.data.sample_data[
            sample_id
        ].hap2level.items():
            if self.taxa_lv in level_dict:
                taxon_name = level_dict[self.taxa_lv]
            elif self.taxa_lv == "haplotype":
                taxon_name = f"{level_dict['species']}_{hap}"
            else:
                raise ValueError(f"Invalid unit_level: {self.taxa_lv}")
            size = self.load_hap_size(sample_id, hap)
            self.taxa_abundance[taxon_name] += (
                size  # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}
            )

    def load_hap_size(self, sample_id: str, hap: str) -> int:
        return int(self.data.sample_data[sample_id].hap_size[hap])

    def normalize_taxa_abundance(self):
        total_size = sum(self.taxa_abundance.values())
        self.taxa_abundance = {
            key: value / total_size * 100
            for key, value in self.taxa_abundance.items()
        }

    def log_taxa_abundance(self):
        self.taxa_abundance = {
            key: np.log(value) for key, value in self.taxa_abundance.items()
        }

    def fill_non_detect_zero(self, taxa_level):
        non_detect_taxa = []
        taxa_set = set(self.taxa_occur[taxa_level])
        sample_id_set = set(self.taxa_occur[SAMPLE_ID_COLUMN])
        for taxa, sample_id in product(taxa_set, sample_id_set):
            if not (
                self.taxa_occur[
                    (self.taxa_occur[taxa_level] == taxa)
                    & (self.taxa_occur[SAMPLE_ID_COLUMN] == sample_id)
                ]
            ).empty:
                continue
            non_detect_taxa.append((taxa, taxa, 0, sample_id))
        self.update_metric_df(
            None,
            non_detect_taxa,
            [
                taxa_level,
                UNIT_COLUMN,
                VALUE_COLUMN,
                SAMPLE_ID_COLUMN,
            ],
            "taxa_occur",
        )

    def convert_taxa_occur_to_taxa_dp(self):
        groupby_columns = self.taxa_occur.columns.drop(
            [UNIT_COLUMN, VALUE_COLUMN, SAMPLE_ID_COLUMN] + self.dp_col
        ).tolist()
        self.dp_df = (
            self.taxa_occur.groupby(groupby_columns)[VALUE_COLUMN]
            .mean()
            .reset_index()
        )
