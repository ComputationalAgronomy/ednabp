import os

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder


class ClusterStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        usearch_prog="usearch",
        in_dir="",
        out_dir="",
        in_suffix=".fasta",
        out_suffix=".fasta",
        minsize: int = 2,
        seq_label: str = "Otu",
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.usearch_prog = usearch_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.cluster_report_suffix = SETTINGS["suffix"]["cluster_report"]
        self.report_suffix = SETTINGS["suffix"]["report"]
        self.parse_params(minsize, seq_label)

    def parse_params(self, minsize, seq_label):
        self.params = f"-minsize {minsize}"
        if seq_label:
            self.params += f" -relabel {seq_label}"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        cluster_outfile = os.path.join(
            self.out_dir, f"{prefix}{self.out_suffix}"
        )
        cluster_report = os.path.join(
            self.out_dir, f"{prefix}{self.cluster_report_suffix}"
        )
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (
            f"{self.usearch_prog} -cluster_otus {self.infile}"
            f" {self.params}"
            f" -otus {cluster_outfile} -uparseout {cluster_report}"
        )
        super().add_stage("Cluster unique sequences into OTUs", cmd)
        super().add_stage_output_to_file(
            "Write usearch report", 0, report, report
        )

    def run(self):
        super().run()
        return all(self.output)


def cluster_demo(config, prefix, derep_dir="", save_dir=""):
    stage = ClusterStage(config, in_dir=derep_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
