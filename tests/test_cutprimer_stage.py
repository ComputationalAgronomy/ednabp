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
        assert stage.in_suffix == ".fastq"
        assert stage.out_suffix == ".fastq"
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
            min_trimmed_read_len=100,
            max_trimmed_read_len=200,
        )

        assert stage.heading == "custom_cutprimer"
        assert stage.cutadapt_prog == "custom_cutadapt"
        assert stage.in_suffix == "_merged.fq"
        assert stage.out_suffix == "_trimmed.fq"

    def test_setup(self, cutprimer_stage, tmp_dirs):
        in_dir, _out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}.fastq")
        with open(input_file, "w") as f:
            f.write("@read1\nACGT\n+\nIIII\n")

        cutprimer_stage.setup(test_prefix)

        assert len(cutprimer_stage.runners) == 2
        command_runner = cutprimer_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}.fastq")
        expected_untrimmed_file = os.path.join(
            _out_dir, "discarded", f"{test_prefix}_untrimmed.fastq"
        )
        # No length limits by default: only empty reads are dropped, and no
        # too-short/too-long files are written
        expected_command = (
            f"cutadapt {expected_infile}"
            f" -g GTCGGTAAAACTCGTGCCAGC;max_error_rate=0.15...CAAACTGGGATTAGATACCCCACTATG;max_error_rate=0.15"
            f" --minimum-length 1 --revcomp -j 4"
            f" --untrimmed-output {expected_untrimmed_file}"
        )

        assert command_runner.command == expected_command
        assert command_runner.prog_name == "Cut primers for merged sequences"


def test_setup_with_length_range(mock_config, tmp_dirs):
    in_dir, out_dir = tmp_dirs
    with open(os.path.join(in_dir, "s1.fastq"), "w") as f:
        f.write("@read1\nACGT\n+\nIIII\n")
    stage = CutPrimerStage(
        config=mock_config,
        in_dir=in_dir,
        out_dir=out_dir,
        min_trimmed_read_len=163,
        max_trimmed_read_len=185,
    )
    stage.setup("s1")

    cmd = stage.runners[0].command
    assert " --minimum-length 163 --maximum-length 185 --revcomp" in cmd
    discarded = os.path.join(out_dir, "discarded", "s1")
    assert cmd.endswith(
        f" --untrimmed-output {discarded}_untrimmed.fastq"
        f" --too-short-output {discarded}_too_short.fastq"
        f" --too-long-output {discarded}_too_long.fastq"
    )
    assert os.path.isdir(os.path.join(out_dir, "discarded"))


def test_setup_with_min_length_only(mock_config, tmp_dirs):
    in_dir, out_dir = tmp_dirs
    with open(os.path.join(in_dir, "s1.fastq"), "w") as f:
        f.write("@read1\nACGT\n+\nIIII\n")
    stage = CutPrimerStage(
        config=mock_config,
        in_dir=in_dir,
        out_dir=out_dir,
        min_trimmed_read_len=100,
    )
    stage.setup("s1")

    cmd = stage.runners[0].command
    assert "--minimum-length 100" in cmd and "--maximum-length" not in cmd
    assert "--too-short-output" in cmd and "--too-long-output" not in cmd
