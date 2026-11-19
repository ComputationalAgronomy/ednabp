import os

import pandas as pd
from Bio import SeqIO

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder


class AddHapStage(stage_builder.StageBuilder):
    """
    Add each OTU/ZOTU's sequence and abundance to the taxonomy table.

    Sequences come from `<prefix><denoise_suffix>` and abundances from the
    OTU table `<prefix><otutab_suffix>` in `denoise_dir`. The table is
    written by the denoise, cluster and cluster_ont stages, so
    `denoise_dir` can be the output directory of any of them.
    """

    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        in_dir="",
        out_dir="",
        in_suffix=".csv",
        out_suffix=".csv",
        denoise_dir="",
        denoise_suffix=".fasta",
        otutab_suffix=None,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.denoise_dir = denoise_dir
        self.denoise_suffix = denoise_suffix
        self.otutab_suffix = otutab_suffix or SETTINGS["suffix"]["otutab"]

    def check_denoise_dir(self):
        if not os.path.isdir(self.denoise_dir):
            raise FileNotFoundError(f"{self.denoise_dir} not found")

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.denoise_file = os.path.join(
            self.denoise_dir, f"{prefix}{self.denoise_suffix}"
        )
        self.otutab_file = os.path.join(
            self.denoise_dir, f"{prefix}{self.otutab_suffix}"
        )
        self.check_infile()
        self.check_denoise_files()
        super().add_stage_function(
            "Add haplotype information", self.add_haplotype_info
        )

    def check_denoise_files(self):
        if not os.path.exists(self.denoise_file):
            raise FileNotFoundError(f"{self.denoise_file} not found")
        if not os.path.exists(self.otutab_file):
            raise FileNotFoundError(f"{self.otutab_file} not found")

    def read_otu_table(self, otutab_path):
        """
        Read a usearch OTU table into {otu: abundance}, summing over all
        sample columns.
        """
        otu_sizes = {}
        with open(otutab_path) as file:
            header = file.readline()
            if not header.startswith("#OTU ID"):
                raise ValueError(
                    f"Unexpected OTU table header in {otutab_path}"
                )
            for line in file:
                fields = line.rstrip("\r\n").split("\t")
                if fields[0]:
                    otu_sizes[fields[0]] = sum(int(x) for x in fields[1:] if x)
        return otu_sizes

    def read_fasta_sequences(self, fasta_path):
        sequences = {}
        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
        return sequences

    def add_haplotype_info(self):
        sequences = {}
        if os.path.exists(self.denoise_file):
            sequences = self.read_fasta_sequences(self.denoise_file)

        zotu_sizes = {}
        if os.path.exists(self.otutab_file):
            zotu_sizes = self.read_otu_table(self.otutab_file)

        all_zotus = set(sequences.keys()) | set(zotu_sizes.keys())

        if os.path.exists(self.infile):
            blast_df = pd.read_csv(self.infile)
        else:
            blast_df = pd.DataFrame()

        zotu_rows = []
        for zotu in all_zotus:
            existing_row = (
                blast_df[blast_df["qseqid"] == zotu]
                if not blast_df.empty
                else pd.DataFrame()
            )
            if existing_row.empty:
                new_row = {"qseqid": zotu}
                zotu_rows.append(new_row)

        if zotu_rows:
            new_df = pd.DataFrame(zotu_rows)
            # Order added rows by the trailing number (Zotu2, Otu10, ...)
            new_df["zotu_num"] = pd.to_numeric(
                new_df["qseqid"].str.extract(r"(\d+)$")[0]
            )
            new_df = new_df.sort_values(by=["zotu_num", "qseqid"])
            new_df = new_df.drop(columns=["zotu_num"])
            blast_df = pd.concat([blast_df, new_df], ignore_index=True)

        blast_df["zotu"] = blast_df["qseqid"].map(sequences).fillna("")
        blast_df["size"] = (
            blast_df["qseqid"].map(zotu_sizes).fillna(0).astype(int)
        )

        blast_df.to_csv(self.outfile, index=False)
        return True

    def run(self):
        super().run()
        return all(self.output)


def addhap_demo(
    config,
    prefix,
    taxa_dir="",
    save_dir="",
    denoise_dir="",
):
    stage = AddHapStage(
        config,
        in_dir=taxa_dir,
        out_dir=save_dir,
        denoise_dir=denoise_dir,
    )
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
