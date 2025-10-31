import os
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from ednabp.bp.step_exec.blast import BlastStage


class TestBlastStage:
    @pytest.fixture
    def blast_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return BlastStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            blast_db="/path/to/blast_db",
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = BlastStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.blast_prog == "blastn"
        assert stage.in_suffix == "_denoise.fasta"
        assert stage.out_suffix == "_blast.csv"

        expected_params = (
            "-db nt -remote -max_target_seqs 1 -evalue 1e-05 "
            "-qcov_hsp_perc 90 -perc_identity 90 "
            '-outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sskingdoms"'
        )
        assert stage.params == expected_params

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = BlastStage(
            config=mock_config,
            heading="custom_blast",
            blast_prog="/path/to/blastn",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_zotus.fa",
            out_suffix="_blast_results.csv",
            blast_db="/path/to/blast_db",
            maxhitnum=5,
            evalue=1e-10,
            qcov_hsp_perc=95,
            perc_identity=97,
            outfmt="6",
            specifiers="qseqid sseqid pident",
        )

        assert stage.heading == "custom_blast"
        assert stage.blast_prog == "/path/to/blastn"
        assert stage.in_suffix == "_zotus.fa"
        assert stage.out_suffix == "_blast_results.csv"

        expected_params = (
            "-db /path/to/blast_db -max_target_seqs 5 -evalue 1e-10 "
            "-qcov_hsp_perc 95 -perc_identity 97 "
            '-outfmt "6 qseqid sseqid pident" -num_threads 4'
        )
        assert stage.params == expected_params

    def test_parse_params_remote_database(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = BlastStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            blast_db="nt",
        )

        assert "nt -remote" in stage.params
        assert "-num_threads" not in stage.params

    def test_setup(self, blast_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_denoise.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n")

        blast_stage.setup(test_prefix)

        assert len(blast_stage.runners) == 2
        blast_runner = blast_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}_denoise.fasta")
        expected_outfile = os.path.join(out_dir, f"{test_prefix}_blast.csv")

        expected_command = (
            f"blastn -query {expected_infile} "
            f"-db /path/to/blast_db -max_target_seqs 1 -evalue 1e-05 "
            f"-qcov_hsp_perc 90 -perc_identity 90 "
            f'-outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sskingdoms" '
            f"-num_threads 4 -out {expected_outfile}"
        )

        assert blast_runner.command == expected_command
        assert blast_runner.prog_name == "Run BLAST"

    @patch("pandas.read_csv")
    def test_add_table_header(self, mock_read_csv, blast_stage, tmp_dirs):
        _in_dir, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        blast_stage.blast_outfile = blast_file

        with open(blast_file, "w") as f:
            f.write("col1,col2,col3\nval1,val2,val3\n")

        mock_blast_df = MagicMock()
        mock_read_csv.return_value = mock_blast_df

        _result = blast_stage.add_table_header()

        mock_read_csv.assert_called_once_with(blast_file)
        mock_blast_df.to_csv.assert_called_once_with(blast_file, index=False)

    @patch("pandas.read_csv")
    def test_add_table_header_empty_result(
        self, mock_read_csv, blast_stage, tmp_dirs
    ):
        _, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        blast_stage.blast_outfile = blast_file

        with open(blast_file, "w") as f:
            f.write("")

        mock_read_csv.side_effect = pd.errors.EmptyDataError("No data")

        result = blast_stage.add_table_header()

        assert result is False
        blast_stage.config.logger.warning.assert_called()