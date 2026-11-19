import os

from ednabp.bp.step_exec.denoise import DenoiseStage
from ednabp.bp.step_exec.utils import read_fasta


def write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def test_per_sample_maps_tagged_copy(mock_config, tmp_dirs):
    """Without a members file, the input is tagged with sample=<prefix>
    and that copy is mapped."""
    in_dir, out_dir = tmp_dirs
    write(os.path.join(in_dir, "s1.fasta"), ">Uniq1;size=10;\nAC\n")
    stage = DenoiseStage(mock_config, in_dir=in_dir, out_dir=out_dir)
    stage.setup("s1")

    query = os.path.join(out_dir, "s1_otutab_query.fasta")
    assert stage.runners[2].prog_name == "Tag sequences with the sample name"
    assert stage.runners[3].command.startswith(f"usearch -otutab {query} ")
    assert stage.runners[4].runner is stage.runners[3]

    stage.runners[2].run()
    assert list(read_fasta(query)) == [("Uniq1;size=10;sample=s1;", "AC")]


def test_pooled_maps_members_file(mock_config, tmp_dirs):
    """After a pool stage, the pool's members file is mapped as is."""
    in_dir, out_dir = tmp_dirs
    write(os.path.join(in_dir, "_pooled.fasta"), ">Uniq1;size=12;\nAC\n")
    members = write(
        os.path.join(in_dir, "_pooled_members.fasta"),
        ">Uniq1;size=10;sample=A;\nAC\n>Uniq1;size=2;sample=B;\nAC\n",
    )
    stage = DenoiseStage(mock_config, in_dir=in_dir, out_dir=out_dir)
    stage.setup("_pooled")

    names = [r.prog_name for r in stage.runners]
    assert "Tag sequences with the sample name" not in names
    assert stage.runners[2].command.startswith(f"usearch -otutab {members} ")
