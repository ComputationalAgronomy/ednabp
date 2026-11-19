import os

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder


class MergeStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        usearch_prog="usearch",
        in_dir="",
        out_dir="",
        in_suffix="_R1.fastq",
        out_suffix=".fastq",
        maxdiff=5,
        pctid=90,
        min_merged_read_len: int | None = None,
        max_merged_read_len: int | None = None,
    ):
        if (
            min_merged_read_len is not None
            and max_merged_read_len is not None
            and min_merged_read_len > max_merged_read_len
        ):
            raise ValueError(
                f"min_merged_read_len ({min_merged_read_len}) is larger than"
                f" max_merged_read_len ({max_merged_read_len})."
            )
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.usearch_prog = usearch_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.report_suffix = SETTINGS["suffix"]["report"]
        self.parse_params(
            maxdiff, pctid, min_merged_read_len, max_merged_read_len
        )

    def parse_params(self, maxdiff, pctid, min_len=None, max_len=None):
        self.params = f"-fastq_maxdiffs {maxdiff} -fastq_pctid {pctid}"
        if min_len is not None:
            self.params += f" -fastq_minmergelen {min_len}"
        if max_len is not None:
            self.params += f" -fastq_maxmergelen {max_len}"
        self.params += f" -threads {self.config.n_cpu}"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        merge_outfile = os.path.join(
            self.out_dir, f"{prefix}{self.out_suffix}"
        )
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (
            f"{self.usearch_prog}"
            f" -fastq_mergepairs {self.infile} -fastqout {merge_outfile}"
            f" {self.params}"
            f" -report {report}"
        )
        super().add_stage("Merge paired-end sequences", cmd)

    def run(self):
        super().run()
        return all(self.output)


def usearch_merge_demo(config, prefix, fastq_dir, save_dir):
    stage = MergeStage(config, in_dir=fastq_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
