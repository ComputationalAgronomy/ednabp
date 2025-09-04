import os

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder


class CutPrimerStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        cutadapt_prog="cutadapt",
        in_dir="",
        out_dir="",
        in_suffix="_merge.fastq",
        out_suffix="_cut.fastq",
        rm_p_5="GTCGGTAAAACTCGTGCCAGC",
        rm_p_3="CAAACTGGGATTAGATACCCCACTATG",
        error_rate=0.15,
        min_read_len=163,
        max_read_len=185,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.cutadapt_prog = cutadapt_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.report_suffix = SETTINGS["suffix"]["report"]
        self.parse_params(
            rm_p_5, rm_p_3, min_read_len, max_read_len, error_rate
        )

    def parse_params(
        self, rm_p_5, rm_p_3, min_read_len, max_read_len, error_rate
    ):
        self.params = (
            f"-g {rm_p_5};max_error_rate={error_rate}...{rm_p_3};max_error_rate={error_rate}"
            f" --minimum-length {min_read_len}"
            f" --maximum-length {max_read_len}"
        )

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        cutprimer_outfile = os.path.join(
            self.out_dir, f"{prefix}{self.out_suffix}"
        )
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (
            f"{self.cutadapt_prog} {self.infile}"
            f" {self.params}"
            f" -j {self.config.n_cpu}"
            f" --untrimmed-output {os.path.join(self.out_dir, prefix + '_untrimmed.fastq')}"
            f" --too-short-output {os.path.join(self.out_dir, prefix + '_too_short.fastq')}"
            f" --too-long-output {os.path.join(self.out_dir, prefix + '_too_long.fastq')}"
        )
        super().add_stage("Cut primers for merged sequences", cmd)
        super().add_stage_output_to_file(
            "Write cutadapt output, report", 0, cutprimer_outfile, report
        )

    def run(self):
        super().run()
        return all(self.output)


def cutadapt_demo(config, prefix, merge_dir="", save_dir=""):
    stage = CutPrimerStage(config, in_dir=merge_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
