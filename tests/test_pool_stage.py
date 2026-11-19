import os

import pytest

from ednabp.bp.step_exec.pool import PoolStage
from ednabp.bp.step_exec.utils import read_fasta


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


@pytest.fixture
def stage(mock_config, tmp_dirs):
    in_dir, out_dir = tmp_dirs
    return PoolStage(mock_config, in_dir=in_dir, out_dir=out_dir)


def test_pool_dereplicated_samples(stage, tmp_dirs):
    """Before OTUs: members carry the size labels; the pooled sequences
    keep size labels (-sizeout); no table is written by the pool."""
    in_dir, out_dir = tmp_dirs
    write(
        os.path.join(in_dir, "A.fasta"),
        ">Uniq1;size=10;\nACGT\n>Uniq2;size=2;\nGGGG\n",
    )
    write(os.path.join(in_dir, "B.fasta"), ">Uniq1;size=30;\nGGGG\n")

    stage.setup(["A", "B"])

    members = os.path.join(out_dir, "_pooled_members.fasta")
    pooled = os.path.join(out_dir, "_pooled.fasta")
    assert stage.pool_otus is False
    assert stage.runners[1].command == (
        f"usearch -fastx_uniques {members} -sizein -sizeout"
        f" -relabel Uniq -fastaout {pooled}"
    )
    assert len(stage.runners) == 3

    stage.write_members()
    assert list(read_fasta(members)) == [
        ("Uniq1;size=10;sample=A;", "ACGT"),
        ("Uniq2;size=2;sample=A;", "GGGG"),
        ("Uniq1;size=30;sample=B;", "GGGG"),
    ]


def test_pool_otus_with_tables(stage, tmp_dirs):
    """After OTUs: counts come from the tables, pooled labels stay plain,
    and the pool writes the table with -otutab -id 1.0."""
    in_dir, out_dir = tmp_dirs
    write(os.path.join(in_dir, "A.fasta"), ">Zotu1\nAAAA\n>Zotu2\nCCCC\n")
    write(os.path.join(in_dir, "A_otu_table.txt"), "#OTU ID\tA\nZotu1\t50\n")
    write(os.path.join(in_dir, "B.fasta"), ">Zotu1\nCCCC\n")
    write(os.path.join(in_dir, "B_otu_table.txt"), "#OTU ID\tB\nZotu1\t70\n")

    stage.setup(["A", "B"])

    members = os.path.join(out_dir, "_pooled_members.fasta")
    pooled = os.path.join(out_dir, "_pooled.fasta")
    assert stage.pool_otus is True
    assert stage.runners[1].command == (
        f"usearch -fastx_uniques {members} -sizein"
        f" -relabel Zotu -fastaout {pooled}"
    )
    assert stage.runners[3].command == (
        f"usearch -otutab {members} -otus {pooled} -id 1.0 -strand plus"
        f" -threads 4 -otutabout "
        + os.path.join(out_dir, "_pooled_otu_table.txt")
    )

    stage.write_members()
    # Zotu2 of A has no count in its table: left out
    assert list(read_fasta(members)) == [
        ("Zotu1;size=50;sample=A;", "AAAA"),
        ("Zotu1;size=70;sample=B;", "CCCC"),
    ]


def test_sample_without_otus_left_out(stage, tmp_dirs):
    """After denoise, a sample with no ZOTUs has an empty FASTA and no OTU
    table: it is left out of the pool with a warning."""
    in_dir, out_dir = tmp_dirs
    write(os.path.join(in_dir, "A.fasta"), ">Zotu1\nAAAA\n")
    write(os.path.join(in_dir, "A_otu_table.txt"), "#OTU ID\tA\nZotu1\t5\n")
    write(os.path.join(in_dir, "B.fasta"), "")

    stage.setup(["A", "B"])

    assert stage.pool_otus is True
    assert stage.prefixes == ["A"]
    warning = stage.config.logger.warning.call_args[0][0]
    assert "B" in warning and "no OTU table" in warning
    stage.write_members()
    members = os.path.join(out_dir, "_pooled_members.fasta")
    assert list(read_fasta(members)) == [("Zotu1;size=5;sample=A;", "AAAA")]


def test_setup_error_is_logged_not_raised(stage):
    """run_all logs a setup error and returns no samples, like other
    stages, instead of stopping the pipeline."""
    assert stage.run_all(["nope"]) == []
    assert (
        "FAIL: setup of stage pool"
        in stage.config.logger.error.call_args[0][0]
    )


def test_custom_pool_name_and_clash(mock_config, tmp_dirs):
    in_dir, out_dir = tmp_dirs
    write(os.path.join(in_dir, "all.fasta"), ">Uniq1;size=1\nAC\n")
    stage = PoolStage(
        mock_config, in_dir=in_dir, out_dir=out_dir, pool_name="all"
    )
    with pytest.raises(ValueError, match="named like the pool"):
        stage.setup(["all"])


def test_missing_input(stage):
    with pytest.raises(FileNotFoundError):
        stage.setup(["nope"])


def test_no_samples(stage):
    assert stage.run_all([]) == []


def test_pipeline_runs_pool_once_for_all_samples(tmp_dirs):
    """The real PoolStage is run once with all samples, and the stages
    after it run on the pool's name."""
    from unittest.mock import patch

    from ednabp.bp.run_bp import BioPipeline

    in_dir, out_dir = tmp_dirs
    for name in ["A", "B"]:
        write(os.path.join(in_dir, f"{name}.fasta"), ">Uniq1;size=3;\nACGT\n")
    with patch.object(BioPipeline, "run_stages_files"):
        pipeline = BioPipeline(
            in_dir, out_dir, enabled_stages=["pool", "denoise"]
        )
    assert PoolStage.pools_samples is True
    pool, denoise = pipeline.stages["pool"], pipeline.stages["denoise"]
    seen = []
    with (
        patch.object(pool, "run_all", return_value=["_pooled"]) as run_all,
        patch.object(
            BioPipeline,
            "run_stage",
            side_effect=lambda name, stage, prefix: seen.append(prefix),
        ),
    ):
        pipeline.run_stages(["A", "B"])
    run_all.assert_called_once_with(["A", "B"])
    assert seen == ["_pooled"]
    assert denoise is pipeline.stages["denoise"]
