import os
import re

import pandas as pd
from Bio import SeqIO

from ..step_build import stage_builder


class AddHapStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        in_dir="",
        out_dir="",
        in_suffix="_taxa.csv",
        out_suffix="_taxa.csv",
        denoise_dir="",
        denoise_suffix="_zotus.fasta",
        report_suffix="_report.txt",
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.denoise_dir = denoise_dir
        self.denoise_suffix = denoise_suffix
        self.report_suffix = re.sub(r"\.\w+", report_suffix, denoise_suffix)

    def check_denoise_dir(self):
        if not os.path.isdir(self.denoise_dir):
            raise FileNotFoundError(f"{self.denoise_dir} not found")

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.denoise_file = os.path.join(
            self.denoise_dir, f"{prefix}{self.denoise_suffix}"
        )
        self.report_file = os.path.join(
            self.denoise_dir, f"{prefix}{self.report_suffix}"
        )
        self.check_infile()
        self.check_denoise_files()
        super().add_stage_function(
            "Add haplotype information", self.add_haplotype_info
        )

    def check_denoise_files(self):
        if not os.path.exists(self.denoise_file):
            raise FileNotFoundError(f"{self.denoise_file} not found")
        if not os.path.exists(self.report_file):
            raise FileNotFoundError(f"{self.report_file} not found")

    def read_denoise_report(self, report_path):
        zotu_sizes = {}
        re_pattern = re.compile(r"size=(\d*)|(zotu)|(chimera)")

        zotus_count = 1
        with open(report_path) as file:
            for line in file:
                if "chfilter" in line and "zotu" in line:
                    matches = ["".join(t) for t in re_pattern.findall(line)]
                    if len(matches) >= 2:
                        size, assigned_type = (
                            matches[0],
                            matches[1],
                        )
                        if assigned_type == "zotu":
                            zotu_sizes[f"Zotu{zotus_count}"] = int(size)
                            zotus_count += 1
        return zotu_sizes

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
        if os.path.exists(self.report_file):
            zotu_sizes = self.read_denoise_report(self.report_file)

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
            new_df["zotu_num"] = (
                new_df["qseqid"].str.extract(r"Zotu(\d+)").astype(int)
            )
            new_df = new_df.sort_values(by="zotu_num")
            new_df = new_df.drop(columns=["zotu_num"])
            blast_df = pd.concat([blast_df, new_df], ignore_index=True)

        blast_df["zotu"] = blast_df["qseqid"].map(sequences).fillna("")
        blast_df["size"] = blast_df["qseqid"].map(zotu_sizes).fillna(0)

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
