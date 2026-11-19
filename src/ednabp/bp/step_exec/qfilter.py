import os

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder


def qscore_to_error_rate(qscore: float) -> float:
    return 10 ** (-qscore / 10)


class QFilterStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        usearch_prog="usearch",
        in_dir="",
        out_dir="",
        in_suffix=".fastq",
        out_suffix=".fasta",
        min_qscore: float = SETTINGS["qfilter"]["min_qscore"],
        max_qscore: int = SETTINGS["qfilter"]["max_qscore"],
        maxns: int | None = SETTINGS["qfilter"]["maxns"],
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.usearch_prog = usearch_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.report_suffix = SETTINGS["suffix"]["report"]
        self.maxee_rate = qscore_to_error_rate(min_qscore)
        self.parse_params(self.maxee_rate, max_qscore, maxns)

    def parse_params(self, maxee_rate, max_qscore, maxns):
        self.params = (
            f"-fastq_maxee_rate {maxee_rate:.6g} -fastq_qmax {max_qscore}"
        )
        if maxns is not None:
            self.params += f" -fastq_maxns {maxns}"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (
            f"{self.usearch_prog} -fastq_filter {self.infile}"
            f" {self.params}"
            f" -fastaout {self.outfile}"
        )
        super().add_stage("Filter low-quality reads", cmd)
        super().add_stage_output_to_file(
            "Write usearch report", 0, report, report
        )

    def run(self):
        super().run()
        return all(self.output)
