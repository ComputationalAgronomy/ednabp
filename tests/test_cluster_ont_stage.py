import os
import random
import shutil

import pytest

from ednabp.bp.step_exec.cluster_ont import (
    ClusterOntStage,
    parse_consensus_labels,
    parse_uc,
)
from ednabp.bp.step_exec.utils import read_fasta


def test_parse_uc(tmp_dirs):
    _, out_dir = tmp_dirs
    uc = os.path.join(out_dir, "c.uc")
    with open(uc, "w") as f:
        f.write(
            "S\t0\t170\t*\t*\t*\t*\t*\tr1\t*\n"
            "H\t0\t170\t98.0\t+\t0\t0\t170M\tr2\tr1\n"
            "S\t1\t170\t*\t*\t*\t*\t*\tr3\t*\n"
            "H\t0\t170\t97.0\t+\t0\t0\t170M\tr4\tr1\n"
            "C\t0\t3\t*\t*\t*\t*\t*\tr1\t*\n"
            "C\t1\t1\t*\t*\t*\t*\t*\tr3\t*\n"
        )
    assert parse_uc(uc) == {0: ["r1", "r2", "r4"], 1: ["r3"]}


def test_parse_consensus_labels(tmp_dirs):
    _, out_dir = tmp_dirs
    fa = os.path.join(out_dir, "cons.fasta")
    with open(fa, "w") as f:
        f.write(">centroid=r1;seqs=3\nACGT\nAC\n>centroid=r3;seqs=1\nTTTT\n")
    assert parse_consensus_labels(fa) == {"r1": "ACGTAC", "r3": "TTTT"}


def test_parse_consensus_labels_annotated(tmp_dirs):
    """Centroid labels with size/sample annotations, as vsearch writes
    them: "centroid=<full label>;seqs=N"."""
    _, out_dir = tmp_dirs
    fa = os.path.join(out_dir, "cons.fasta")
    with open(fa, "w") as f:
        f.write(">centroid=S2:Uniq1;size=7;sample=S2;;seqs=2\nACGT\n")
    assert parse_consensus_labels(fa) == {"S2:Uniq1;size=7;sample=S2;": "ACGT"}


def test_read_fasta_multiline(tmp_dirs):
    _, out_dir = tmp_dirs
    fa = os.path.join(out_dir, "x.fasta")
    with open(fa, "w") as f:
        f.write(">a comment\nAC\nGT\n\n>b\nTT\n")
    assert list(read_fasta(fa)) == [("a comment", "ACGT"), ("b", "TT")]


class TestClusterOntStage:
    def test_init_defaults(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = ClusterOntStage(mock_config, in_dir=in_dir, out_dir=out_dir)
        assert stage.vsearch_prog == "vsearch"
        assert stage.cluster_id == 0.90
        assert stage.min_cluster_size == 10
        assert stage.racon_rounds == 2
        assert stage.medaka is False
        assert stage.consensus_id == 0.99
        assert stage.otutab_id == 0.90
        assert stage.cluster_strand is False
        assert stage.otutab_strand is True

    def test_setup(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        infile = os.path.join(in_dir, "s1.fasta")
        with open(infile, "w") as f:
            f.write(">r1\nACGT\n")
        stage = ClusterOntStage(mock_config, in_dir=in_dir, out_dir=out_dir)
        stage.setup("s1")

        tmp = os.path.join(out_dir, "s1_cluster_ont_tmp")
        assert stage.runners[0].command == (
            f"vsearch --cluster_fast {infile} --id 0.9 --strand plus"
            f" --threads 4 --uc {tmp}/clusters.uc"
            f" --consout {tmp}/draft_consensus.fasta"
        )
        assert stage.runners[2].prog_name == (
            "Polish cluster consensus sequences"
        )
        assert stage.runners[3].prog_name == (
            "Tag sequences with the sample name"
        )
        query = os.path.join(out_dir, "s1_otutab_query.fasta")
        otutab_cmd = stage.runners[4].command
        assert otutab_cmd.startswith(f"usearch -otutab {query} -otus ")
        assert "-otus " + os.path.join(out_dir, "s1.fasta") in otutab_cmd
        # clustering defaults to plus strand, the OTU table to both
        assert "-id 0.9 -strand both" in otutab_cmd
        # reads that match no OTU go to the discarded subfolder
        notmatched = os.path.join(
            out_dir, "discarded", "s1_otu_notmatched.fasta"
        )
        assert otutab_cmd.endswith(f" -notmatched {notmatched}")
        assert os.path.isdir(os.path.join(out_dir, "discarded"))
        # vsearch + report, polish, tag + otutab + report, tmp cleanup
        assert len(stage.runners) == 7

    def test_setup_separate_strands(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        with open(os.path.join(in_dir, "s1.fasta"), "w") as f:
            f.write(">r1\nACGT\n")
        stage = ClusterOntStage(
            mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            cluster_strand=True,
            otutab_strand=False,
        )
        stage.setup("s1")
        assert "--strand both" in stage.runners[0].command
        assert "-strand plus" in stage.runners[4].command

    def test_setup_keep_tmp(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        with open(os.path.join(in_dir, "s1.fasta"), "w") as f:
            f.write(">r1\nACGT\n")
        stage = ClusterOntStage(
            mock_config, in_dir=in_dir, out_dir=out_dir, keep_tmp=True
        )
        stage.setup("s1")
        assert len(stage.runners) == 6

    def test_setup_missing_infile(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = ClusterOntStage(mock_config, in_dir=in_dir, out_dir=out_dir)
        with pytest.raises(FileNotFoundError):
            stage.setup("missing")


def simulate_reads(path, seed=1):
    """Two templates with ~3% ONT-like errors (mostly indels)."""
    rng = random.Random(seed)
    templates = {
        name: "".join(rng.choice("ACGT") for _ in range(170))
        for name in ("A", "B")
    }

    def mutate(s):
        out = []
        for c in s:
            r = rng.random()
            if r < 0.012:
                continue
            if r < 0.024:
                out += [c, c]
                continue
            if r < 0.03:
                c = rng.choice([x for x in "ACGT" if x != c])
            out.append(c)
        return "".join(out)

    with open(path, "w") as f:
        n = 0
        for name, count in (("A", 40), ("B", 20)):
            for _ in range(count):
                n += 1
                f.write(f">r{n}_{name} runid=x\n{mutate(templates[name])}\n")
    return templates


@pytest.mark.skipif(
    not all(
        shutil.which(p) for p in ("vsearch", "minimap2", "racon", "usearch")
    ),
    reason="vsearch, minimap2, racon and usearch must be on PATH",
)
def test_cluster_ont_end_to_end(tmp_dirs):
    from unittest.mock import Mock

    from ednabp.common.config import Config

    in_dir, out_dir = tmp_dirs
    templates = simulate_reads(os.path.join(in_dir, "s1.fasta"))
    config = Mock(spec=Config)
    config.verbose, config.dry, config.n_cpu = False, False, 2
    config.logger = Mock()

    stage = ClusterOntStage(
        config, in_dir=in_dir, out_dir=out_dir, min_cluster_size=5
    )
    stage.setup("s1")
    assert stage.run()

    otus = dict(read_fasta(os.path.join(out_dir, "s1.fasta")))
    assert sorted(otus.values()) == sorted(templates.values())
    assert not os.path.exists(os.path.join(out_dir, "s1_cluster_ont_tmp"))
    with open(os.path.join(out_dir, "s1_otu_table.txt")) as f:
        rows = dict(line.split() for line in f.readlines()[1:])
    assert int(rows["Otu1"]) >= 38 and int(rows["Otu2"]) >= 18
