from abc import ABC
from typing import TYPE_CHECKING

from . import base_logger

if TYPE_CHECKING:
    from ..sample_data import SampleData


class Writer(ABC):
    def __init__(self, sampledata: "SampleData", verbose: bool):
        self.logger = base_logger.logger
        if verbose:
            self.logger.setLevel("INFO")
        else:
            self.logger.setLevel("WARNING")
        self.sample_id_used = None
        self._import_data(sampledata)

    def _import_data(self, samplesdata):
        self.sample_data = samplesdata.sample_data
        self.sample_id_list = samplesdata.sample_id_list
        if samplesdata.sample_metadata:
            self.sample_metadata = samplesdata.sample_metadata
        if samplesdata.spc_info:
            self.spc_info = samplesdata.spc_info

    @base_logger.prog_log(prog_name="Load sample ID list")
    def _load_sample_id_list(self, sample_id_list: list[str] | None = None):
        if sample_id_list is None:
            self.logger.info(
                f"No sample ID list specified. Using all {len(self.sample_id_list)} samples."
            )
            self.sample_id_used = self.sample_id_list
        else:
            self.logger.info(f"Specified {len(sample_id_list)} samples.")
            for sample_id in sample_id_list:
                if sample_id not in self.sample_id_list:
                    raise ValueError(
                        f"Specified invalid sample ID: {sample_id}."
                    )
            self.sample_id_used = sample_id_list

    def _add_file_handler(self, log_path):
        fh = base_logger.get_file_handler(log_path)
        self.logger.addHandler(fh)
