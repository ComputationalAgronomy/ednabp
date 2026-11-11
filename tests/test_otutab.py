import os

import pytest

from ednabp.bp.step_exec.otutab import sum_row


def test_collapse_otutab(tmp_dirs):
    _, out_dir = tmp_dirs
    path = os.path.join(out_dir, "s1_otu_table.txt")
    with open(path, "w") as f:
        f.write(
            "#OTU ID\tUniq3\tUniq1\tUniq2\nOtu2\t0\t0\t6\nOtu1\t3\t10\t0\n"
        )

    sum_row(path, "sample001")

    with open(path) as f:
        assert f.read() == "#OTU ID\tsample001\nOtu2\t6\nOtu1\t13\n"


def test_collapse_otutab_single_column(tmp_dirs):
    _, out_dir = tmp_dirs
    path = os.path.join(out_dir, "t.txt")
    with open(path, "w") as f:
        f.write("#OTU ID\tUniq1\nZotu1\t40\n")

    sum_row(path, "s1")

    with open(path) as f:
        assert f.read() == "#OTU ID\ts1\nZotu1\t40\n"


def test_collapse_otutab_bad_header(tmp_dirs):
    _, out_dir = tmp_dirs
    path = os.path.join(out_dir, "t.txt")
    with open(path, "w") as f:
        f.write("not a table\n")

    with pytest.raises(ValueError):
        sum_row(path, "s1")
