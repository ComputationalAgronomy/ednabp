from abc import ABC

from . import base_logger
from ..runner_exec import data_container


class Writer(ABC):
    def __init__(self, sampledata: data_container.SampleData, no_verbose: bool):
        self.logger = base_logger.logger
        if no_verbose:
            self.logger.setLevel("WARNING")
        self.sample_id_used = None
        self._import_data(sampledata)

    def _import_data(self, samplesdata):
        self.sample_data = samplesdata.sample_data
        self.sample_id_list = samplesdata.sample_id_list
        if hasattr(samplesdata, "sample_metadata"):
            self.sample_metadata = samplesdata.sample_metadata

    @base_logger.prog_log(prog_name="Load sample ID list")
    def _load_sample_id_list(self, sample_id_list: list[str] | None = None):
        if sample_id_list is None:
            self.logger.info(f"No sample ID list specified. Using all {len(self.sample_id_list)} samples.")
            self.sample_id_used = self.sample_id_list
        else:
            self.logger.info(f"Specified {len(sample_id_list)} samples.")
            for sample_id in sample_id_list:
                if sample_id not in self.sample_id_list:
                    raise ValueError(f"Specified invalid sample ID: {sample_id}.")
            self.sample_id_used = sample_id_list

    def _add_file_handler(self, log_path):
        fh = base_logger.get_file_handler(log_path)
        self.logger.addHandler(fh)





