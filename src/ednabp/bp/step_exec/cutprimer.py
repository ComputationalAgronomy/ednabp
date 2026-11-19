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
        in_suffix=".fastq",
        out_suffix=".fastq",
        rm_p_5="GTCGGTAAAACTCGTGCCAGC",
        rm_p_3="CAAACTGGGATTAGATACCCCACTATG",
        error_rate=0.15,
        min_trimmed_read_len: int | None = None,
        max_trimmed_read_len: int | None = None,
        revcomp: bool = True,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.cutadapt_prog = cutadapt_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.report_suffix = SETTINGS["suffix"]["report"]
        self.min_len = min_trimmed_read_len
        self.max_len = max_trimmed_read_len
        self.parse_params(
            rm_p_5,
            rm_p_3,
            min_trimmed_read_len,
            max_trimmed_read_len,
            error_rate,
            revcomp,
        )

    def parse_params(
        self, rm_p_5, rm_p_3, min_len, max_len, error_rate, revcomp
    ):
        self.params = f"-g {rm_p_5};max_error_rate={error_rate}...{rm_p_3};max_error_rate={error_rate}"
        self.params += f" --minimum-length {min_len if min_len else 1}"
        if max_len is not None:
            self.params += f" --maximum-length {max_len}"
        if revcomp:
            self.params += " --revcomp"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        cutprimer_outfile = os.path.join(
            self.out_dir, f"{prefix}{self.out_suffix}"
        )
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        discarded_dir = os.path.join(self.out_dir, SETTINGS["discarded_dir"])
        os.makedirs(discarded_dir, exist_ok=True)
        discarded = os.path.join(discarded_dir, prefix)
        ext = self.out_suffix
        cmd = (
            f"{self.cutadapt_prog} {self.infile}"
            f" {self.params}"
            f" -j {self.config.n_cpu}"
            f" --untrimmed-output {discarded}_untrimmed{ext}"
        )
        if self.min_len is not None:
            cmd += f" --too-short-output {discarded}_too_short{ext}"
        if self.max_len is not None:
            cmd += f" --too-long-output {discarded}_too_long{ext}"
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
