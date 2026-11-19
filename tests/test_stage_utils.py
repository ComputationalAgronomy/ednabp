import os

import pytest

from ednabp.bp.step_exec.utils import (
    label_size,
    members,
    read_table,
    seq_id,
)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


@pytest.mark.parametrize(
    "label, sid, size",
    [
        ("Uniq1;size=10;", "Uniq1", 10),
        ("Uniq1;size=10", "Uniq1", 10),
        ("Zotu3", "Zotu3", 1),
        ("r1 runid=x ch=5", "r1", 1),
        ("r1;size=4; runid=x", "r1", 4),
    ],
)
def test_labels(label, sid, size):
    assert seq_id(label) == sid
    assert label_size(label) == size


def test_read_table(tmp_dirs):
    _, d = tmp_dirs
    path = write(
        os.path.join(d, "t.txt"), "#OTU ID\tA\tB\nOtu1\t3\t0\nOtu2\t0\t7\n"
    )
    assert read_table(path) == (["A", "B"], {"Otu1": [3, 0], "Otu2": [0, 7]})


def test_read_table_bad_header(tmp_dirs):
    _, d = tmp_dirs
    path = write(os.path.join(d, "t.txt"), "Uniq1;size=87;\tchfilter zotu\n")
    with pytest.raises(ValueError):
        read_table(path)


def test_members_from_size_labels(tmp_dirs):
    _, d = tmp_dirs
    fasta = write(
        os.path.join(d, "s.fa"),
        ">Uniq1;size=10;\nAC\n>r2 runid=x\nAG\n",
    )
    assert list(members(fasta, "S01")) == [
        ("Uniq1;size=10;sample=S01;", "AC"),
        ("r2;size=1;sample=S01;", "AG"),
    ]


def test_members_from_counts_skip_zero(tmp_dirs):
    _, d = tmp_dirs
    fasta = write(os.path.join(d, "z.fa"), ">Zotu1\nAC\n>Zotu2\nAG\n")
    assert list(members(fasta, "A", {"Zotu1": 7})) == [
        ("Zotu1;size=7;sample=A;", "AC")
    ]


@pytest.mark.parametrize("sample", ["", "a b", "a;b", "a=b"])
def test_members_bad_sample_name(tmp_dirs, sample):
    _, d = tmp_dirs
    fasta = write(os.path.join(d, "z.fa"), ">Zotu1\nAC\n")
    with pytest.raises(ValueError):
        list(members(fasta, sample))
