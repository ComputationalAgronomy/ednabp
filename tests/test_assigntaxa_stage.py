import csv
import os
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from ednabp.bp.step_exec.assigntaxa import AssignTaxaStage, get_empty_rank


class TestAssignTaxaStage:
    @pytest.fixture
    def sample_lineage_data(self):
        return [
            {
                "genus_name": "Salmo",
                "family_name": "Salmonidae",
                "order_name": "Salmoniformes",
                "class_name": "Actinopterygii",
                "phylum_name": "Chordata",
                "kingdom_name": "Animalia",
            },
            {
                "genus_name": "Oncorhynchus",
                "family_name": "Salmonidae",
                "order_name": "Salmoniformes",
                "class_name": "Actinopterygii",
                "phylum_name": "Chordata",
                "kingdom_name": "Animalia",
            },
        ]

    @pytest.fixture
    def lineage_file(self, tmp_dirs, sample_lineage_data):
        _in_dir, out_dir = tmp_dirs
        lineage_path = os.path.join(out_dir, "lineage.csv")

        with open(lineage_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=sample_lineage_data[0].keys()
            )
            writer.writeheader()
            writer.writerows(sample_lineage_data)

        return lineage_path

    @pytest.fixture
    def assigntaxa_stage(self, mock_config, tmp_dirs, lineage_file):
        in_dir, out_dir = tmp_dirs
        return AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            blast_db="/path/to/blast_db",
            lineage_db=lineage_file,
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
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
        assert stage.lineage_db == "nucleotide"

        expected_params = (
            "-db nt -remote -max_target_seqs 1 -evalue 1e-05 "
            "-qcov_hsp_perc 90 -perc_identity 90 "
            '-outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sskingdoms"'
        )
        assert stage.params == expected_params

    def test_init_custom(self, mock_config, tmp_dirs, lineage_file):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config,
            heading="custom_blast",
            blast_prog="/path/to/blastn",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_zotus.fa",
            out_suffix="_taxonomy.csv",
            blast_db="/path/to/blast_db",
            lineage_db=lineage_file,
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
        assert stage.out_suffix == "_taxonomy.csv"
        assert stage.lineage_db == lineage_file

        expected_params = (
            "-db /path/to/blast_db -max_target_seqs 5 -evalue 1e-10 "
            "-qcov_hsp_perc 95 -perc_identity 97 "
            '-outfmt "6 qseqid sseqid pident" -num_threads 4'
        )
        assert stage.params == expected_params

    def test_parse_params_remote_database(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            blast_db="nt",
            lineage_db=None,
        )

        assert "nt -remote" in stage.params
        assert "-num_threads" not in stage.params

    def test_parse_genus2otherlv(self, assigntaxa_stage, sample_lineage_data):
        expected_mapping = {
            "Salmo": [
                "Animalia",
                "Chordata",
                "Actinopterygii",
                "Salmoniformes",
                "Salmonidae",
            ],
            "Oncorhynchus": [
                "Animalia",
                "Chordata",
                "Actinopterygii",
                "Salmoniformes",
                "Salmonidae",
            ],
        }

        assert hasattr(assigntaxa_stage, "genus2otherlv")
        assert assigntaxa_stage.genus2otherlv == expected_mapping

    def test_parse_lineage_db_invalid_file(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with pytest.raises(ValueError):
            AssignTaxaStage(
                config=mock_config,
                in_dir=in_dir,
                out_dir=out_dir,
                blast_db="/path/to/db",
                lineage_db="/nonexistent/path.csv",
            )

    def test_setup(self, assigntaxa_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_denoise.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n")

        assigntaxa_stage.setup(test_prefix)

        assert len(assigntaxa_stage.runners) >= 2
        blast_runner = assigntaxa_stage.runners[0]

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
        assert len(assigntaxa_stage.runners) == 3

    def test_setup_no_lineage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        assigntaxa_stage = AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            lineage_db=None,
        )

        input_file = os.path.join(in_dir, "test_denoise.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGT\n")

        assigntaxa_stage.setup("test")

        assert len(assigntaxa_stage.runners) == 2

    @patch("pandas.read_csv")
    def test_add_table_header(self, mock_read_csv, assigntaxa_stage, tmp_dirs):
        _in_dir, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assigntaxa_stage.blast_outfile = blast_file

        with open(blast_file, "w") as f:
            f.write("col1,col2,col3\nval1,val2,val3\n")

        mock_blast_df = MagicMock()
        mock_read_csv.return_value = mock_blast_df

        _result = assigntaxa_stage.add_table_header()

        mock_read_csv.assert_called_once_with(blast_file)
        mock_blast_df.to_csv.assert_called_once_with(blast_file, index=False)

    @patch("pandas.read_csv")
    def test_add_table_header_empty_result(
        self, mock_read_csv, assigntaxa_stage, tmp_dirs
    ):
        _, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assigntaxa_stage.blast_outfile = blast_file

        with open(blast_file, "w") as f:
            f.write("")

        mock_read_csv.side_effect = pd.errors.EmptyDataError("No data")

        result = assigntaxa_stage.add_table_header()

        assert result is False
        assigntaxa_stage.config.logger.warning.assert_called()

    @patch("pandas.read_csv")
    def test_get_lineage_from_custom_db(
        self, mock_read_csv, assigntaxa_stage, tmp_dirs
    ):
        _, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assigntaxa_stage.blast_outfile = blast_file

        mock_blast_df = pd.DataFrame(
            {"sseqid": ["seq1|Salmo_trutta", "seq2|Unknown_species"]}
        )
        mock_read_csv.return_value = mock_blast_df

        with patch.object(mock_blast_df, "merge") as mock_merge:
            mock_merged = MagicMock()
            mock_merge.return_value = mock_merged

            result = assigntaxa_stage.get_lineage_from_custom_db()

            assert result is True
            mock_merged.to_csv.assert_called_once_with(blast_file, index=False)

    def test_get_empty_rank(self):
        result = get_empty_rank()
        expected = {
            "kingdom": "",
            "phylum": "",
            "class": "",
            "order": "",
            "family": "",
            "genus": "",
            "species": "",
        }
        assert result == expected

    def test_nucleotide_lineage_db(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            blast_db="/path/to/db",
            lineage_db="nucleotide",
        )

        assert stage.lineage_db == "nucleotide"
        assert hasattr(stage, "get_lineage_func")
        assert stage.entrez_email == "your.email@example.com"
