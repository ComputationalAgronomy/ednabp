import os
import tempfile
from unittest.mock import Mock, patch

import pytest

from ednabp.bp.step_exec.dereplicate import DereplicateStage
from ednabp.common.config import Config


class TestDereplicateStage:
    @pytest.fixture
    def dereplicate_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return DereplicateStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = DereplicateStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.usearch_prog == "usearch"
        assert stage.in_suffix == "_cut.fasta"
        assert stage.out_suffix == "_uniq.fasta"
        assert stage.report_suffix == "_report.txt"
        assert stage.write_report is True

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = DereplicateStage(
            config=mock_config,
            heading="custom_derep",
            usearch_prog="usearch11",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_trimmed.fa",
            out_suffix="_unique.fa",
            annot_size=False,
            seq_label="Custom",
            write_report=False,
        )

        assert stage.heading == "custom_derep"
        assert stage.usearch_prog == "usearch11"
        assert stage.in_suffix == "_trimmed.fa"
        assert stage.out_suffix == "_unique.fa"
        assert stage.write_report is False

    def test_setup(self, dereplicate_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_cut.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGT\n>seq2\nACGT\n>seq3\nTGCA\n")

        dereplicate_stage.setup(test_prefix)

        assert len(dereplicate_stage.runners) == 2
        command_runner = dereplicate_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}_cut.fasta")
        expected_outfile = os.path.join(out_dir, f"{test_prefix}_uniq.fasta")

        expected_command = (
            f"usearch -fastx_uniques {expected_infile} "
            f"-sizeout -relabel Uniq -threads 4 "
            f"-fastaout {expected_outfile}"
        )

        assert command_runner.command == expected_command
        assert command_runner.prog_name == "Dereplicate trimmed sequences"
