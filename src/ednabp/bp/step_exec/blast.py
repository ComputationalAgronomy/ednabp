import os
from typing import Literal

import pandas as pd

from ..step_build import stage_builder


class BlastStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        blast_prog="blastn",
        in_dir="",
        out_dir="",
        in_suffix="_denoise.fasta",
        out_suffix="_blast.csv",
        maxhitnum: int = 1,
        evalue: float = 0.00001,
        qcov_hsp_perc: int = 90,
        perc_identity: int = 90,
        outfmt: str = "10",
        specifiers: str = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        blast_db: str | Literal["nt"] = "nt",
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.blast_prog = blast_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.blast_outfile = None
        self.specifiers = specifiers
        self.parse_params(
            maxhitnum,
            evalue,
            qcov_hsp_perc,
            perc_identity,
            outfmt,
            blast_db,
        )

    def parse_params(
        self,
        maxhitnum,
        evalue,
        qcov_hsp_perc,
        perc_identity,
        outfmt,
        blast_db,
    ):
        if blast_db == "nt":
            blast_db = "nt -remote"
        self.params = (
            f"-max_target_seqs {maxhitnum}"
            f" -evalue {evalue}"
            f" -qcov_hsp_perc {qcov_hsp_perc}"
            f" -perc_identity {perc_identity}"
            f' -outfmt "{outfmt} {self.specifiers}"'
            f" -db {blast_db}"
        )
        if "remote" not in self.params:
            self.params += f" -num_threads {self.config.n_cpu}"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.blast_outfile = os.path.join(
            self.out_dir, f"{prefix}{self.out_suffix}"
        )
        self.check_infile()
        cmd = (
            f"{self.blast_prog} -query {self.infile}"
            f" {self.params}"
            f" -out {self.blast_outfile}"
        )
        super().add_stage("Run BLAST", cmd)
        super().add_stage_function(
            "Update output table header", self.add_table_header
        )

    def add_table_header(self):
        if os.path.exists(self.blast_outfile):
            try:
                blast_result = pd.read_csv(self.blast_outfile, header=None)
                blast_result.columns = self.specifiers.split(" ")
                blast_result.to_csv(self.blast_outfile, index=False)
            except pd.errors.EmptyDataError:
                self.config.logger.warning(
                    f"BLAST result is empty: {self.blast_outfile}"
                )
                return False
        else:
            self.config.logger.error(
                f"BLAST result does not exist: {self.blast_outfile}"
            )
            return False

    def run(self):
        super().run()
        return all(self.output)


def blast_demo(
    config,
    prefix,
    denoise_dir="",
    save_dir="",
    db_path="",
    specifiers="qseqid sscinames sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
):
    stage = BlastStage(
        config,
        in_dir=denoise_dir,
        out_dir=save_dir,
        blast_db=db_path,
        specifiers=specifiers,
    )
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
