from Bio import SeqIO
import os
from typing import override

from fastq_processor.step_build import stage_builder


class FqToFaStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_fq_to_fa.py", cutprimer_dir="", save_dir="",
                 overwrite=True
        ):
        super().__init__(heading=heading, config=config)
        self.in_suffix = "cut.fastq"
        self.out_suffix = "cut.fasta"
        self.cutprimer_dir = cutprimer_dir
        self.save_dir = save_dir

    def fq_to_fa(self):
        with open(self.infile, "r") as fq, open(self.outfile, "w") as fa:
            for record in SeqIO.parse(fq, "fastq"):
                SeqIO.write(record, fa, "fasta")

    def setup(self, prefix):
        self.infile = os.path.join(self.cutprimer_dir, f"{prefix}_{self.in_suffix}")
        self.outfile = os.path.join(self.save_dir, f"{prefix}_{self.out_suffix}")
        self.check_infile()
        self.check_savedir()
        super().add_stage_function("Convert FASTQ to FASTA", self.fq_to_fa)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, cutprimer_dir="", save_dir=""):
    stage = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete