import os
import re

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder
from .utils import members, read_fasta, read_table, seq_id, write_fasta


def label_prefix(fasta_path):
    if not os.path.exists(fasta_path):  # dry run
        return "Seq"
    for label, _ in read_fasta(fasta_path):
        return re.match(r"\D*", seq_id(label)).group() or "Seq"
    return "Seq"


class PoolStage(stage_builder.StageBuilder):
    pools_samples = True

    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        usearch_prog="usearch",
        in_dir="",
        out_dir="",
        in_suffix=".fasta",
        out_suffix=".fasta",
        pool_name: str = "_pooled",
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.usearch_prog = usearch_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.pool_name = pool_name
        self.suffix = SETTINGS["suffix"]

    def in_path(self, prefix, suffix):
        return os.path.join(self.in_dir, f"{prefix}{suffix}")

    def out_path(self, suffix):
        return os.path.join(self.out_dir, f"{self.pool_name}{suffix}")

    def setup(self, prefixes):
        if isinstance(prefixes, str):
            prefixes = [prefixes]
        self.prefixes = list(prefixes)
        if self.pool_name in self.prefixes:
            raise ValueError(
                f"A sample is named like the pool ({self.pool_name});"
                " set another pool_name."
            )
        for prefix in self.prefixes:
            self.infile = self.in_path(prefix, self.in_suffix)
            self.check_infile()

        no_table = [
            p
            for p in self.prefixes
            if not os.path.exists(self.in_path(p, self.suffix["otutab"]))
        ]
        self.pool_otus = len(no_table) < len(self.prefixes)
        if self.pool_otus and no_table:
            self.config.logger.warning(
                f"Not pooled, no OTU table (no OTUs): {', '.join(no_table)}"
            )
            self.prefixes = [p for p in self.prefixes if p not in no_table]

        self.outfile = self.out_path(self.out_suffix)
        self.members_file = self.out_path(self.suffix["members"])
        self.otutab_file = self.out_path(self.suffix["otutab"])
        report = self.out_path(self.suffix["report"])
        otutab_report = self.out_path(self.suffix["otutab_report"])
        relabel = label_prefix(self.in_path(self.prefixes[0], self.in_suffix))
        self.config.logger.info(
            f"Pooling {'OTUs' if self.pool_otus else 'sequences'} of"
            f" {len(self.prefixes)} samples"
        )

        super().add_stage_function(
            f"Tag and combine {len(self.prefixes)} samples",
            self.write_members,
        )
        sizeout = "" if self.pool_otus else " -sizeout"
        super().add_stage(
            "Merge identical sequences",
            f"{self.usearch_prog} -fastx_uniques {self.members_file}"
            f" -sizein{sizeout} -relabel {relabel}"
            f" -fastaout {self.outfile}",
        )
        super().add_stage_output_to_file(
            "Write usearch report", len(self.runners) - 1, report, report
        )
        if self.pool_otus:
            super().add_stage(
                "Map sample OTUs to pooled OTUs",
                f"{self.usearch_prog} -otutab {self.members_file}"
                f" -otus {self.outfile} -id 1.0 -strand plus"
                f" -threads {self.config.n_cpu}"
                f" -otutabout {self.otutab_file}",
            )
            super().add_stage_output_to_file(
                "Write otutab report",
                len(self.runners) - 1,
                otutab_report,
                otutab_report,
            )

    def run_all(self, prefixes):
        if not prefixes:
            return []
        try:
            self.setup(prefixes)
        except (FileNotFoundError, ValueError) as e:
            self.config.logger.error(f"FAIL: setup of stage pool: {e}")
            self.runners = []
            return []
        return [self.pool_name] if self.run() else []

    def sample_members(self, prefix):
        fasta = self.in_path(prefix, self.in_suffix)
        sizes = None
        if self.pool_otus:
            _, rows = read_table(self.in_path(prefix, self.suffix["otutab"]))
            sizes = {sid: sum(counts) for sid, counts in rows.items()}
        return members(fasta, prefix, sizes)

    def write_members(self):
        write_fasta(
            (r for p in self.prefixes for r in self.sample_members(p)),
            self.members_file,
        )

    def run(self):
        super().run()
        return all(self.output)
