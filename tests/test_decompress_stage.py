import gzip
import os

import pytest

from ednabp.bp.step_exec.decompress import DecompressStage


class TestDecompressStage:
    @pytest.fixture
    def decompress_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return DecompressStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

    @pytest.fixture
    def sample_fastq_content(self):
        return (
            b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nIIIIIIII\n"
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = DecompressStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.insuffix_list == ["_R1.fastq.gz", "_R2.fastq.gz"]
        assert stage.outsuffix_list == ["_R1.fastq", "_R2.fastq"]
        assert stage.infile_list == []
        assert stage.outfile_list == []

    def test_setup_existed_files(
        self, decompress_stage, tmp_dirs, sample_fastq_content
    ):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        for suffix in ["_R1.fastq.gz", "_R2.fastq.gz"]:
            filepath = os.path.join(in_dir, f"{test_prefix}{suffix}")
            with gzip.open(filepath, "wb") as f:
                f.write(sample_fastq_content)

        decompress_stage.setup(test_prefix)

        assert len(decompress_stage.infile_list) == 2
        assert len(decompress_stage.outfile_list) == 2

        expected_infiles = [
            os.path.join(in_dir, f"{test_prefix}_R1.fastq.gz"),
            os.path.join(in_dir, f"{test_prefix}_R2.fastq.gz"),
        ]
        expected_outfiles = [
            os.path.join(out_dir, f"{test_prefix}_R1.fastq"),
            os.path.join(out_dir, f"{test_prefix}_R2.fastq"),
        ]

        assert decompress_stage.infile_list == expected_infiles
        assert decompress_stage.outfile_list == expected_outfiles

        assert len(decompress_stage.runners) == 1
        runner = decompress_stage.runners[0]
        assert runner.prog_name == "Decompress FASTQ.GZ to FASTQ"

    def test_setup_missing_file(self, decompress_stage):
        with pytest.raises(
            FileNotFoundError, match="nonexistent_R1.fastq.gz not found"
        ):
            decompress_stage.setup("nonexistent")

    def test_run(self, decompress_stage, tmp_dirs, sample_fastq_content):
        in_dir, _out_dir = tmp_dirs
        test_prefix = "test_sample"

        input_files = []
        for suffix in ["_R1.fastq.gz", "_R2.fastq.gz"]:
            filepath = os.path.join(in_dir, f"{test_prefix}{suffix}")
            with gzip.open(filepath, "wb") as f:
                f.write(sample_fastq_content)
            input_files.append(filepath)

        decompress_stage.setup(test_prefix)

        decompress_stage.run()

        for outfile in decompress_stage.outfile_list:
            assert os.path.exists(outfile)
            with open(outfile, "rb") as f:
                content = f.read()
                assert content == sample_fastq_content

    def test_run_empty_file(self, decompress_stage, tmp_dirs):
        in_dir, _out_dir = tmp_dirs
        test_prefix = "empty_sample"

        for suffix in ["_R1.fastq.gz", "_R2.fastq.gz"]:
            filepath = os.path.join(in_dir, f"{test_prefix}{suffix}")
            with gzip.open(filepath, "wb") as f:
                f.write(b"")

        decompress_stage.setup(test_prefix)
        decompress_stage.run()

        for outfile in decompress_stage.outfile_list:
            assert os.path.exists(outfile)
            with open(outfile, "rb") as f:
                assert f.read() == b""
