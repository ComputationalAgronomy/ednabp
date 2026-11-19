import os

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder
from .utils import members, write_fasta


def add_otutab_stages(
    stage: stage_builder.StageBuilder,
    prefix: str,
    uniq_file: str,
    otu_file: str,
    zotus: bool = False,
    identity: float = 0.97,
    strand: bool = True,
):
    suffix = SETTINGS["suffix"]
    out = {
        key: os.path.join(stage.out_dir, f"{prefix}{suffix[key]}")
        for key in [
            "otutab",
            "otu_mapping",
            "otu_notmatched",
            "otutab_report",
            "otutab_query",
        ]
    }

    discarded_dir = os.path.join(stage.out_dir, SETTINGS["discarded_dir"])
    os.makedirs(discarded_dir, exist_ok=True)
    out["otu_notmatched"] = os.path.join(
        discarded_dir, f"{prefix}{suffix['otu_notmatched']}"
    )

    query = os.path.join(stage.in_dir, f"{prefix}{suffix['members']}")
    if not os.path.exists(query):
        query = out["otutab_query"]
        stage.add_stage_function(
            "Tag sequences with the sample name",
            lambda: write_fasta(members(uniq_file, prefix), query),
        )

    db_option = "-zotus" if zotus else "-otus"
    strand_option = "both" if strand else "plus"
    cmd = (
        f"{stage.usearch_prog} -otutab {query}"
        f" {db_option} {otu_file}"
        f" -id {identity} -strand {strand_option}"
        f" -threads {stage.config.n_cpu}"
        f" -otutabout {out['otutab']} -mapout {out['otu_mapping']}"
        f" -notmatched {out['otu_notmatched']}"
    )
    stage.add_stage("Map sequences to OTUs", cmd)
    stage.add_stage_output_to_file(
        "Write otutab report",
        len(stage.runners) - 1,
        out["otutab_report"],
        out["otutab_report"],
    )
