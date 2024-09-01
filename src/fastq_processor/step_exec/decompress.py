import os
import gzip
from typing import override

from fastq_processor.step_build import stage_builder


class DecompressStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_gzip_decompress.py",
                 in_dir="", out_dir="",
                 in_suffix="R1.fastq.gz", out_suffix="R1.fastq"
        ):
        super().__init__(heading=heading, config=config, in_dir=in_dir, out_dir=out_dir)
        in_suffix2 = in_suffix.replace("R1", "R2")
        out_suffix2 = out_suffix.replace("R1", "R2")
        self.insuffix_list = [in_suffix, in_suffix2]
        self.outsuffix_list = [out_suffix, out_suffix2]
        self.infile_list = []
        self.outfile_list = []

    def gunzip(self):
        for infile, outfile in zip(self.infile_list, self.outfile_list):
            with gzip.open(infile, 'rb')as in_handler, open(outfile, 'wb') as out_handler:
                for line in in_handler:
                    out_handler.write(line)

    def setup(self, prefix):
        for in_suffix in self.insuffix_list:
            self.infile = os.path.join(self.in_dir, f"{prefix}_{in_suffix}")
            self.check_infile()
            self.infile_list.append(self.infile)

        for out_suffix in self.outsuffix_list:
            self.outfile = os.path.join(self.out_dir, f"{prefix}_{out_suffix}")
            self.outfile_list.append(self.outfile)
        super().add_stage_function("Decompress FASTQ.GZ to FASTQ", self.gunzip)

    def run(self):
        super().run()
        return all(self.output)


def decompress_demo(config, prefix, fastq_dir="", save_dir=""):
    stage = DecompressStage(config, fastq_dir, save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete