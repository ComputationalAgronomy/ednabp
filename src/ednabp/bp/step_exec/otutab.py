import os

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder


def sum_row(otutab_file, sample):
    rows = []
    with open(otutab_file) as f:
        header = f.readline()
        if not header.startswith("#OTU ID"):
            raise ValueError(f"Unexpected OTU table header in {otutab_file}")
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if fields[0]:
                rows.append((fields[0], sum(int(x) for x in fields[1:] if x)))

    with open(otutab_file, "w") as f:
        f.write(f"#OTU ID\t{sample}\n")
        for otu, count in rows:
            f.write(f"{otu}\t{count}\n")


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
    otutab_file = os.path.join(stage.out_dir, f"{prefix}{suffix['otutab']}")
    mapping_file = os.path.join(
        stage.out_dir, f"{prefix}{suffix['otu_mapping']}"
    )
    notmatched_file = os.path.join(
        stage.out_dir, f"{prefix}{suffix['otu_notmatched']}"
    )
    otutab_report = os.path.join(
        stage.out_dir, f"{prefix}{suffix['otutab_report']}"
    )

    db_option = "-zotus" if zotus else "-otus"
    strand_option = "both" if strand else "plus"
    cmd = (
        f"{stage.usearch_prog} -otutab {uniq_file}"
        f" {db_option} {otu_file}"
        f" -id {identity} -strand {strand_option}"
        f" -threads {stage.config.n_cpu}"
        f" -otutabout {otutab_file} -mapout {mapping_file}"
        f" -notmatched {notmatched_file}"
    )
    stage.add_stage("Generate OTU abundance table", cmd)
    stage.add_stage_output_to_file(
        "Write otutab report",
        len(stage.runners) - 1,
        otutab_report,
        otutab_report,
    )
    stage.add_stage_function(
        "Sum OTU table rows into one sample column",
        lambda: sum_row(otutab_file, prefix),
    )
