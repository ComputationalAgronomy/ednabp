from abc import ABC

from ..common import base_logger, config
from ..data import BPData, MitoData


class BaseWriter(ABC):
    def __init__(self, data: BPData | MitoData, verbose: bool):
        self.sample_id_used = None
        self.config = config.Config()
        self.config.verbose = verbose
        self.import_data(data)

    def import_data(self, data):
        if not isinstance(data, BPData, MitoData):
            raise ValueError("data must be an instance of BPData or MitoData")
        if not hasattr(data, "sample_data"):
            raise ValueError(
                "No attribute `sample_data` found in the input data."
            )
        if not hasattr(data, "sample_id_list"):
            raise ValueError(
                "No attribute `sample_id_list` found in the input data."
            )
        if not hasattr(data, "sample_metadata"):
            self.config.logger.warning(
                "WARNING: No attribute `sample_metadata` found in the input data."
            )
        if not hasattr(data, "spc_info"):
            self.config.logger.warning(
                "WARNING: No attribute `spc_info` found in the input data."
            )
        self.data = data

    @base_logger.prog_log(prog_name="Load sample ID list")
    def load_sample_id_list(self, sample_id_list: list[str] | None = None):
        if sample_id_list is None:
            self.config.logger.info(
                f"No sample ID list specified. Using all {len(self.data.sample_id_list)} samples."
            )
            self.sample_id_used = self.data.sample_id_list
        else:
            self.config.logger.info(
                f"Specified {len(sample_id_list)} samples."
            )
            for sample_id in sample_id_list:
                if sample_id not in self.data.sample_id_list:
                    raise ValueError(
                        f"Specified invalid sample ID: {sample_id}."
                    )
            self.sample_id_used = sample_id_list
