import os

from fastq_processor.step_build import stage_builder


class DenoiseStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_usearch_denoise.py",
                 in_dir="", out_dir="",
                 in_suffix="_uniq.fasta", out_suffix="_denoise.fasta",
                 minsize=8,
                 alpha=2
                 ):
        super().__init__(heading=heading, config=config, in_dir=in_dir, out_dir=out_dir)
        self.USEARCH_PROG = "usearch" # TODO(SW): Don't use `.exe`, doesn't make sense in docker/ubuntu
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.denoise_report_suffix = "_denoise_report.txt"
        self.report_suffix = "_report.txt"
        self.parse_params(minsize, alpha)

    def parse_params(self, minsize, alpha):
        self.params = (f"-minsize {minsize} -unoise_alpha {alpha}")

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        denoise_outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        denoise_report = os.path.join(self.out_dir, f"{prefix}{self.denoise_report_suffix}")
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (f"{self.USEARCH_PROG} -unoise3 {self.infile}"
               f" {self.params}"
               f" -threads {self.config.n_cpu}"
               f" -zotus {denoise_outfile} -tabbedout {denoise_report}")
        super().add_stage("Denoise unique sequences", cmd)
        super().add_stage_output_to_file("Write usearch report", 0, report, report)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, derep_dir="", save_dir=""):
    stage = DenoiseStage(config, in_dir=derep_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete