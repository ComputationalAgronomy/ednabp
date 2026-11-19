import os

import pandas as pd
from Bio import SeqIO

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder
from .utils import read_table


class AddHapStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        in_dir="",
        out_dir="",
        in_suffix=".csv",
        out_suffix=".csv",
        otu_dir="",
        otu_suffix=".fasta",
        otutab_suffix=None,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.otu_dir = otu_dir
        self.otu_suffix = otu_suffix
        self.otutab_suffix = otutab_suffix or SETTINGS["suffix"]["otutab"]

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.otu_file = os.path.join(
            self.otu_dir, f"{prefix}{self.otu_suffix}"
        )
        self.otutab_file = os.path.join(
            self.otu_dir, f"{prefix}{self.otutab_suffix}"
        )
        self.check_infile()
        self.check_otu_files()
        super().add_stage_function(
            "Add haplotype information", self.add_haplotype_info
        )

    def check_otu_files(self):
        if self.config.dry:
            return
        if not os.path.exists(self.otu_file):
            raise FileNotFoundError(f"{self.otu_file} not found")
        if not os.path.exists(self.otutab_file):
            raise FileNotFoundError(f"{self.otutab_file} not found")

    def read_otu_table(self, otutab_path):
        samples, rows = read_table(otutab_path)
        return samples, {
            otu: dict(zip(samples, counts, strict=True))
            for otu, counts in rows.items()
        }

    def add_sample_columns(self, blast_df, samples, otu_counts):
        clash = sorted(set(samples) & set(blast_df.columns))
        if clash:
            raise ValueError(
                f"Sample names clash with table columns: {', '.join(clash)}"
            )
        blast_df = blast_df.drop(columns=["size"], errors="ignore")
        counts = pd.DataFrame(
            [
                [otu_counts.get(otu, {}).get(s, 0) for s in samples]
                for otu in blast_df["qseqid"]
            ],
            columns=samples,
            index=blast_df.index,
            dtype=int,
        )
        return pd.concat([blast_df, counts], axis=1)

    def read_fasta_sequences(self, fasta_path):
        sequences = {}
        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
        return sequences

    def add_haplotype_info(self):
        otu_sequences = {}
        if os.path.exists(self.otu_file):
            otu_sequences = self.read_fasta_sequences(self.otu_file)

        samples, otu_counts = [], {}
        if os.path.exists(self.otutab_file):
            samples, otu_counts = self.read_otu_table(self.otutab_file)

        if os.path.exists(self.infile):
            blast_df = pd.read_csv(self.infile)
        else:
            blast_df = pd.DataFrame()

        otu_ids = set(otu_sequences) | set(otu_counts)
        with_hit = set(blast_df["qseqid"]) if not blast_df.empty else set()
        missing_rows = [{"qseqid": otu_id} for otu_id in otu_ids - with_hit]
        if missing_rows:
            new_df = pd.DataFrame(missing_rows)
            new_df["otu_number"] = pd.to_numeric(
                new_df["qseqid"].str.extract(r"(\d+)$")[0]
            )
            new_df = new_df.sort_values(by=["otu_number", "qseqid"])
            new_df = new_df.drop(columns=["otu_number"])
            blast_df = pd.concat([blast_df, new_df], ignore_index=True)

        blast_df = blast_df.drop(columns=["sequence"], errors="ignore")
        blast_df["sequence"] = blast_df["qseqid"].map(otu_sequences).fillna("")
        blast_df = self.add_sample_columns(blast_df, samples, otu_counts)

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
    otu_dir="",
):
    stage = AddHapStage(
        config,
        in_dir=taxa_dir,
        out_dir=save_dir,
        otu_dir=otu_dir,
    )
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
