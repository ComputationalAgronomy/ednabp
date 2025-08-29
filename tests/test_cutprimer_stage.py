import os

import pytest

from ednabp.bp.step_exec.cutprimer import CutPrimerStage


class TestCutPrimerStage:
    @pytest.fixture
    def cutprimer_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return CutPrimerStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = CutPrimerStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.cutadapt_prog == "cutadapt"
        assert stage.in_suffix == "_merge.fastq"
        assert stage.out_suffix == "_cut.fastq"
        assert stage.report_suffix == "_report.txt"

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = CutPrimerStage(
            config=mock_config,
            heading="custom_cutprimer",
            cutadapt_prog="custom_cutadapt",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_merged.fq",
            out_suffix="_trimmed.fq",
            rm_p_5="AAAA",
            rm_p_3="TTTT",
            error_rate=0.2,
            min_read_len=100,
            max_read_len=200,
        )

        assert stage.heading == "custom_cutprimer"
        assert stage.cutadapt_prog == "custom_cutadapt"
        assert stage.in_suffix == "_merged.fq"
        assert stage.out_suffix == "_trimmed.fq"

    def test_setup(self, cutprimer_stage, tmp_dirs):
        in_dir, _out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_merge.fastq")
        with open(input_file, "w") as f:
            f.write("@read1\nACGT\n+\nIIII\n")

        cutprimer_stage.setup(test_prefix)

        assert len(cutprimer_stage.runners) == 2
        command_runner = cutprimer_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}_merge.fastq")
        expected_command = (
            f"cutadapt {expected_infile} "
            f"-g GTCGGTAAAACTCGTGCCAGC;max_error_rate=0.15...CAAACTGGGATTAGATACCCCACTATG;max_error_rate=0.15 "
            f"--minimum-length 163 --maximum-length 185 --discard-untrimmed -j 4"
        )

        assert command_runner.command == expected_command
        assert command_runner.prog_name == "Cut primers for merged sequences"
