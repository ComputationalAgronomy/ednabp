import os

from fastq_processor.step_build import stage_builder


class DereplicateStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_usearch_dereplicate.py",
                 in_dir="", out_dir="",
                 in_suffix="_cut.fasta", out_suffix="_uniq.fasta",
                 annot_size: bool = True,
                 seq_label: str = "Uniq"
                 ):
        super().__init__(heading=heading, config=config, in_dir=in_dir, out_dir=out_dir)
        self.USEARCH_PROG = "usearch" # TODO(SW): Don't use `.exe`, doesn't make sense in docker/ubuntu
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.report_suffix = "_report.txt"
        self.parse_params(annot_size, seq_label)

    def parse_params(self, annot_size, seq_label):
        sizeout = "-sizeout" if annot_size else ""
        self.params = (f"{sizeout} -relabel {seq_label}")

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        dereplicate_outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (f"{self.USEARCH_PROG} -fastx_uniques {self.infile}"
               f" {self.params}"
               f" -threads {self.config.n_cpu}"
               f" -fastaout {dereplicate_outfile}")
        super().add_stage("Dereplicate trimmed sequences", cmd)
        super().add_stage_output_to_file("Write usearch report", 0, report, report)

    def run(self):
        super().run()
        return all(self.output)


def dereplicate_demo(config, prefix, fasta_dir="", save_dir=""):
    stage = DereplicateStage(config, in_dir=fasta_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete