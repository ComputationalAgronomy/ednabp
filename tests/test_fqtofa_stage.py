import os

import pytest

from ednabp.bp.step_exec.fqtofa import FqToFaStage


class TestFqToFaStage:
    @pytest.fixture
    def fqtofa_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return FqToFaStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

    @pytest.fixture
    def sample_fastq_content(self):
        return "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nIIIIIIII\n"

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = FqToFaStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.in_suffix == ".fastq"
        assert stage.out_suffix == ".fasta"

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = FqToFaStage(
            config=mock_config,
            heading="custom_fqtofa",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_trimmed.fq",
            out_suffix="_trimmed.fa",
        )

        assert stage.heading == "custom_fqtofa"
        assert stage.in_suffix == "_trimmed.fq"
        assert stage.out_suffix == "_trimmed.fa"

    def test_setup(self, fqtofa_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}.fastq")
        with open(input_file, "w") as f:
            f.write("@read1\nACGT\n+\nIIII\n")

        fqtofa_stage.setup(test_prefix)

        assert len(fqtofa_stage.runners) == 1
        runner = fqtofa_stage.runners[0]
        assert runner.prog_name == "Convert FASTQ to FASTA"

        expected_infile = os.path.join(in_dir, f"{test_prefix}.fastq")
        expected_outfile = os.path.join(out_dir, f"{test_prefix}.fasta")
        assert fqtofa_stage.infile == expected_infile
        assert fqtofa_stage.outfile == expected_outfile

    def test_run(self, fqtofa_stage, tmp_dirs, sample_fastq_content):
        in_dir, out_dir = tmp_dirs
        test_prefix = "real_test"

        input_file = os.path.join(in_dir, f"{test_prefix}.fastq")
        with open(input_file, "w") as f:
            f.write(sample_fastq_content)

        fqtofa_stage.setup(test_prefix)

        fqtofa_stage.run()

        output_file = os.path.join(out_dir, f"{test_prefix}.fasta")
        assert os.path.exists(output_file)

        with open(output_file) as f:
            content = f.read()
            assert content.startswith(">")
            assert "ACGTACGT" in content
            assert "TGCATGCA" in content
            assert "@" not in content
            assert "+" not in content

    def test_empty_file(self, fqtofa_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "empty_test"

        input_file = os.path.join(in_dir, f"{test_prefix}.fastq")
        with open(input_file, "w") as f:
            f.write("")

        fqtofa_stage.setup(test_prefix)
        fqtofa_stage.fq_to_fa()

        output_file = os.path.join(out_dir, f"{test_prefix}.fasta")
        assert os.path.exists(output_file)

        with open(output_file) as f:
            content = f.read()
            assert content == ""
