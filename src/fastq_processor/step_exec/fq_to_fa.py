from Bio import SeqIO
import os

from fastq_processor.step_build import stage_builder


class FqToFaStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_fq_to_fa.py",
                 in_dir="", out_dir="",
                 in_suffix="cut.fastq", out_suffix="cut.fasta"
        ):
        super().__init__(heading=heading, config=config, in_dir=in_dir, out_dir=out_dir)
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix

    def fq_to_fa(self):
        with open(self.infile, "r") as fq, open(self.outfile, "w") as fa:
            for record in SeqIO.parse(fq, "fastq"):
                SeqIO.write(record, fa, "fasta")

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}_{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}_{self.out_suffix}")
        self.check_infile()
        super().add_stage_function("Convert FASTQ to FASTA", self.fq_to_fa)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, cutprimer_dir="", save_dir=""):
    stage = FqToFaStage(config, in_dir=cutprimer_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete