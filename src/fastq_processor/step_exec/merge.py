import os

from fastq_processor.step_build import stage_builder


class MergeStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_usearch_merge.py",
                 in_dir="", out_dir="",
                 in_suffix="_R1.fastq", out_suffix="_merge.fastq",
                 maxdiff=5,
                 pctid=90
                 ):
        super().__init__(heading=heading, config=config, in_dir=in_dir, out_dir=out_dir)
        self.USEARCH_PROG = "usearch" # TODO(SW): Don't use `.exe`, doesn't make sense in docker/ubuntu
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.report_suffix = "_report.txt"
        self.parse_params(maxdiff, pctid)

    def parse_params(self, maxdiff, pctid):
        self.params = (f"-fastq_maxdiffs {maxdiff} -fastq_pctid {pctid}"
                       f" -threads {self.config.n_cpu}")

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        merge_outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        cmd = (f"{self.USEARCH_PROG}"
               f" -fastq_mergepairs {self.infile} -fastqout {merge_outfile}"
               f" {self.params}"
               f" -report {report}")
        super().add_stage("Merge paired-end sequences", cmd)

    def run(self):
        super().run()
        return all(self.output)


def usearch_merge_demo(config, prefix, fastq_dir, save_dir):
    stage = MergeStage(config, in_dir=fastq_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete