import os

import pytest

from ednabp.bp.step_exec.denoise import DenoiseStage


class TestDenoiseStage:
    @pytest.fixture
    def denoise_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return DenoiseStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = DenoiseStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.usearch_prog == "usearch"
        assert stage.in_suffix == ".fasta"
        assert stage.out_suffix == ".fasta"
        assert stage.denoise_report_suffix == "_denoise_report.txt"
        assert stage.report_suffix == "_report.txt"

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = DenoiseStage(
            config=mock_config,
            heading="custom_denoise",
            usearch_prog="usearch12",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_unique.fa",
            out_suffix="_zotus.fa",
            minsize=10,
            alpha=3,
        )

        assert stage.heading == "custom_denoise"
        assert stage.usearch_prog == "usearch12"
        assert stage.in_suffix == "_unique.fa"
        assert stage.out_suffix == "_zotus.fa"

    def test_setup(self, denoise_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGT\n>seq2\nTGCA\n")

        denoise_stage.setup(test_prefix)

        # Command + output redirection, tag + otutab + report
        assert len(denoise_stage.runners) == 5
        command_runner = denoise_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}.fasta")
        expected_outfile = os.path.join(out_dir, f"{test_prefix}.fasta")
        expected_denoise_report = os.path.join(
            out_dir, f"{test_prefix}_denoise_report.txt"
        )

        expected_command = (
            f"usearch -unoise3 {expected_infile} "
            f"-minsize 8 -unoise_alpha 2 -threads 4 "
            f"-zotus {expected_outfile} -tabbedout {expected_denoise_report}"
        )

        assert command_runner.command == expected_command
        assert command_runner.prog_name == "Denoise unique sequences"

        otutab_runner = denoise_stage.runners[3]
        expected_otutab_command = (
            f"usearch -otutab "
            f"{os.path.join(out_dir, f'{test_prefix}_otutab_query.fasta')} "
            f"-zotus {expected_outfile} -id 0.97 -strand both -threads 4 "
            f"-otutabout {os.path.join(out_dir, f'{test_prefix}_otu_table.txt')} "
            f"-mapout {os.path.join(out_dir, f'{test_prefix}_otu_mapping.txt')} "
            f"-notmatched {os.path.join(out_dir, 'discarded', f'{test_prefix}_otu_notmatched.fasta')}"
        )
        assert otutab_runner.command == expected_otutab_command
