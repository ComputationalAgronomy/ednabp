import os

import pytest

from ednabp.bp.step_exec.merge import MergeStage


class TestMergeStage:
    @pytest.fixture
    def merge_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return MergeStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

    def test_stage_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = MergeStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.USEARCH_PROG == "usearch"
        assert stage.in_suffix == "_R1.fastq"
        assert stage.out_suffix == "_merge.fastq"
        assert stage.report_suffix == "_report.txt"
        assert stage.params == "-fastq_maxdiffs 5 -fastq_pctid 90 -threads 4"

    def test_stage_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = MergeStage(
            config=mock_config,
            heading="custom_merge",
            usearch_prog="usearch12",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_R1.fq",
            out_suffix="_merged.fq",
            maxdiff=10,
            pctid=85,
        )

        assert stage.heading == "custom_merge"
        assert stage.USEARCH_PROG == "usearch12"
        assert stage.in_suffix == "_R1.fq"
        assert stage.out_suffix == "_merged.fq"
        assert stage.params == "-fastq_maxdiffs 10 -fastq_pctid 85 -threads 4"

    def test_setup_creates_correct_command(self, merge_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_R1.fastq")
        with open(input_file, "w") as f:
            f.write("@test\nACGT\n+\nIIII\n")

        merge_stage.setup(test_prefix)

        assert len(merge_stage.runners) == 1
        runner = merge_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}_R1.fastq")
        expected_outfile = os.path.join(out_dir, f"{test_prefix}_merge.fastq")
        expected_report = os.path.join(out_dir, f"{test_prefix}_report.txt")

        expected_command = (
            f"usearch -fastq_mergepairs {expected_infile} "
            f"-fastqout {expected_outfile} "
            f"-fastq_maxdiffs 5 -fastq_pctid 90 -threads 4 "
            f"-report {expected_report}"
        )

        assert runner.command == expected_command
        assert runner.prog_name == "Merge paired-end sequences"
