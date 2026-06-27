import os

import pandas as pd

from ..step_build import stage_builder
from .addlineage import RANKS


def normalize_species_names(name):
    if pd.isna(name) or name == "":
        return []
    if "_x_" in name:
        return name.split("_x_")
    parts = name.split("_")
    if len(parts) >= 3:
        return [name, "_".join(parts[:2])]
    if "_sp." in name:
        return [name, name.split("_")[0] + "_sp."]
    return [name]


def filter_by_score(group, score_column, tol_pct):
    if group.empty:
        return group
    highest = group[score_column].max()
    return group[group[score_column] >= highest * (1 - tol_pct / 100)]


def lca_of_lineages(lineages):
    result = [""] * len(RANKS)
    lca_rank = "NA"
    last_i = len(RANKS) - 1

    for i, col in enumerate(RANKS):
        per_hit = []
        any_missing = False
        for lineage in lineages:
            candidates = set(normalize_species_names(lineage[i]))
            if candidates:
                per_hit.append(candidates)
            else:
                any_missing = True

        if any_missing:
            if i == last_i:
                break
            continue

        if not per_hit:
            break

        common = per_hit[0].copy()
        for s in per_hit[1:]:
            common &= s

        if common:
            result[i] = max(common, key=len)
            lca_rank = col
        else:
            break

    return result, lca_rank


class LcaStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        in_dir="",
        out_dir="",
        in_suffix=".csv",
        out_suffix=".csv",
        tol_pct: float = 1.0,
        score_column: str = "bitscore",
        qseqid_column: str = "qseqid",
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.tol_pct = tol_pct
        self.score_column = score_column
        self.qseqid_column = qseqid_column

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.check_infile()
        if not self.config.dry:
            super().add_stage_function("Run LCA", self.run_lca)

    def run_lca(self):
        try:
            df = pd.read_csv(self.infile)
            non_tax_cols = [c for c in df.columns if c not in RANKS]
            results = []

            for seq_name, group in df.groupby(self.qseqid_column, sort=False):
                filtered = filter_by_score(
                    group, self.score_column, self.tol_pct
                )
                if filtered.empty:
                    continue

                lineages = filtered[RANKS].values.tolist()
                lca_lineage, lca_rank = lca_of_lineages(lineages)

                best_row = filtered.loc[filtered[self.score_column].idxmax()]
                row = {col: best_row[col] for col in non_tax_cols}
                for i, col in enumerate(RANKS):
                    row[col] = lca_lineage[i]
                row["lca_rank"] = lca_rank
                row["num_assignments"] = len(filtered)
                results.append(row)

            pd.DataFrame(results).to_csv(self.outfile, index=False)
            return True
        except Exception as e:
            self.config.logger.error(f"LCA failed: {e}")
            return False

    def run(self):
        super().run()
        return all(self.output)
