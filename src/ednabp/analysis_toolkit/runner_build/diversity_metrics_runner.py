from abc import ABC
from collections import defaultdict
import os

import pandas as pd

# from . import base_runner, base_logger
from ednabp.analysis_toolkit.runner_build import base_runner, base_logger


class DMRunner(base_runner.Runner):
    '''
    An abstract class for running diversity metrics related analysis
    Metrics include:
        - abundance: sum of sequence read for each target taxon. (e.g. How many sequences of fish got in sampleA?)
        - richness: number of unique units for each target taxon. (e.g. How many different fish species got in sampleA?)
    '''
    ABUNDANCE_COLUMN = "Abundance"
    RICHNESS_COLUMN = "Richness"
    UNIT_COLUMN = "Unit"
    SAMPLE_ID_COLUMN = "Sample_id"

    def __init__(self, sampledata, no_verbose):
        super().__init__(sampledata, no_verbose)

    @base_logger.prog_log("Write richness table")
    def write_richness_table(self,
            save_dir: str,
            taxa_level: str,
            unit_level: str = "species",
            sample_id_list: list[str] | None = None,
        ):
        """
        Write the richness data to a CSV file.

        :param save_dir: The directory path to save the CSV file.
        :param taxa_level: The taxa level to group units (e.g., species, family, etc.).
        :param unit_level: The taxa level to calculate richness (e.g., species, haplotype, etc.). Default is "species".
        :param sample_id_list: List of sample IDs to import. Default: None (all available sample IDs will be used).
        """

        self._load_sample_id_list(sample_id_list)
        self._create_richness_df(taxa_level, unit_level)
        DMRunner._export_df(save_dir, f"{taxa_level}_richness.csv", self.richness_df)
        self.analysis_type = "Write species richness to csv"

    @base_logger.prog_log("Write Abundance table")
    def write_abundance_table(self,
            save_dir: str,
            taxa_level: str,
            normalize: bool = False,
            sample_id_list: list[str] | None = None
        ):
        """
        Write the abundance data to a CSV file.

        :param save_dir: The directory to save the CSV file.
        :param taxa_level: The name of the level to target (e.g., species, family, etc.).
        :param normalize: Whether to normalize the abundance data.
        :param sample_id_list: A list of sample IDs to write.
        """
        os.makedirs(save_dir, exist_ok=True)

        self._load_sample_id_list(sample_id_list)
        self._create_abundance_df(taxa_level, normalize)
        DMRunner._export_df(save_dir, f"{taxa_level}_abundance.csv", self.abundance_df)
        self.analysis_type = "Write species abundance to csv"

    @base_logger.prog_log("Calculate taxa richness and create dataframe")
    def _create_richness_df(self, taxa_level, unit_level):
        self.occurrence_df = pd.DataFrame()
        for sample_id in self.sample_id_used:
            self.logger.info(f"Sample ID: {sample_id}")
            self._get_sample_units_occurence(sample_id, taxa_level, unit_level)
            self._update_metric_df(sample_id,
                                   self.units_occurrence,
                                   [taxa_level, self.UNIT_COLUMN, self.RICHNESS_COLUMN],
                                   "occurrence_df") 
        # self._filter_richness_df()
        self._convert_units_occurrence_to_taxa_richness(taxa_level) # columns: taxa_level, Richness, Sample_id
        self._add_sample_metadata("richness_df") # columns: taxa_level, Richness, Sample_id, Site, Year, Month, Sample

    @base_logger.prog_log("Calculate taxa abundance and create dataframe")
    def _create_abundance_df(self, taxa_level: str, normalize: bool):
        self.abundance_df = pd.DataFrame()
        for sample_id in self.sample_id_used:
            self.logger.info(f"Sample ID: {sample_id}")
            self._get_sample_taxa_abundance(sample_id, taxa_level)
            if normalize:
                self._normalize_taxa_abundance()
            self._update_metric_df(sample_id,
                                   self.taxa_abundance.items(),
                                   [taxa_level, self.ABUNDANCE_COLUMN],
                                   "abundance_df") # columns: taxa_level, Abundance, Sample_id
        self._add_sample_metadata("abundance_df") # columns: taxa_level, Abundance, Sample_id, Site, Year, Month, Sample

    @base_logger.prog_log("Export dataframe to CSV file")
    def _export_df(save_dir, file_name, metric_df):
        os.makedirs(save_dir, exist_ok=True)
        output_path = os.path.join(save_dir, file_name)
        metric_df.to_csv(output_path, index=False)

    def _get_sample_units_occurence(self, sample_id: str, taxa_level: str, unit_level: str):
        """
        Get the appear units for a given sample.
        """
        self.units_occurrence = []
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            assert taxa_level in level_dict, f"Invalid taxa_level: {taxa_level}"
            taxon_name = level_dict[taxa_level]
            unit_name = DMRunner._get_unit_name(level_dict, unit_level, hap)
            unit = (taxon_name, unit_name, 1)
            if unit not in self.units_occurrence:
                self.units_occurrence.append(unit) # e.g. {'SpA': 1, 'SpB': 1, 'SpC': 1}

    def _get_unit_name(level_dict: dict, unit_level: str, hap: str) -> str:
        if unit_level in level_dict:
            return level_dict[unit_level]
        if unit_level == "haplotype":
            return f"{level_dict['species']}_{hap}"
        raise ValueError(f"Invalid unit_level: {unit_level}")

    def _update_metric_df(self, *args):
        sample_id, metric_data, metric_columns_name, metric_df_name = args
        df = pd.DataFrame(metric_data, columns=metric_columns_name)
        df[self.SAMPLE_ID_COLUMN] = sample_id
        metric_df = getattr(self, metric_df_name)
        updated_metric_df = pd.concat([metric_df, df], ignore_index=True)
        setattr(self, metric_df_name, updated_metric_df)
    
    def _add_sample_metadata(self, metric_df_name):
        metric_df = getattr(self, metric_df_name)
        metadata_df = (pd.DataFrame
            .from_dict(self.sample_metadata, orient="index")
            .reset_index()
            .rename(columns={'index': self.SAMPLE_ID_COLUMN})
        )
        updated_metric_df = pd.merge(
            left=metric_df,
            right=metadata_df,
            on=self.SAMPLE_ID_COLUMN,
            how="outer"
        )
        if updated_metric_df.isna().any():
            self.logger.warning(f"WARNING: Some samples are missing metadata. Filling them with 'Unknown'.")
            updated_metric_df = updated_metric_df.fillna("Unknown")
        setattr(self, metric_df_name, updated_metric_df)

    # def _filter_richness_df(self, site_occur_thres: int = 0, sample_occur_thres: int = 0):
    #     sp_list = self.richness_df.groupby(["Unit"])[["Site","Year","Month","Sample"]].nunique()
    #     sp_list = sp_list[sp_list["Site"]>0] # Filtering by number of sites
    #     sp_list = sp_list[sp_list["Sample"]>0] # Filtering by number of samples
    #     survive = np.array([])
    #     for m in sp_list.index:
    #         survive = np.append(survive, np.where(self.richness_df[sp_list.index.name]==m)[0])
    #     self.richness_df = self.richness_df.iloc[survive,:].reset_index(drop=True)

    def _convert_units_occurrence_to_taxa_richness(self, taxa_level):
        self.richness_df = (
            self.occurrence_df.groupby([taxa_level, self.SAMPLE_ID_COLUMN])
            [self.RICHNESS_COLUMN]
            .sum()
            .reset_index()
        )

    def _get_sample_taxa_abundance(self, sample_id: str, taxa_level: str):
        self.taxa_abundance = defaultdict(int)
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            if taxa_level in level_dict:
                taxon_name = level_dict[taxa_level]
            elif taxa_level == "haplotype":
                taxon_name = f"{level_dict["species"]}_{hap}"
            else:
                raise ValueError(f"Invalid unit_level: {taxa_level}")
            size = self._load_hap_size(sample_id, hap)
            self.taxa_abundance[taxon_name] += size # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}

    def _load_hap_size(self, sample_id: str, hap: str) -> int:
        return int(self.sample_data[sample_id].hap_size[hap])

    def _normalize_taxa_abundance(self):
        total_size = sum(self.taxa_abundance.values())
        self.taxa_abundance = {key: value/total_size * 100 for key, value in self.taxa_abundance.items()}

    def _save_html(self, fig_type: str, save_html_dir: str, save_name: str):
        """
        Save the barchart as an HTML file.

        :param save_html_dir: The directory to save the HTML file.
        :param save_html_name: The name of the HTML file. If not provided, the name will be "{level}_barchart". Default is None.
        """
        fig_path = os.path.join(save_html_dir, f"{save_name}.html")
        self.fig.write_html(fig_path)
        self.logger.info(f"{fig_type} saved to: {fig_path}")