import csv
import os
from unittest.mock import MagicMock, Mock, patch

import pandas as pd
import pytest

from ednabp.bp.step_exec.assigntaxa import AssignTaxaStage


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
    def assign_taxa_stage(self, mock_config, tmp_dirs, lineage_file):
        in_dir, out_dir = tmp_dirs
        return AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            db_path="/path/to/database",
            lineage_path=lineage_file,
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            db_path="/path/to/db",
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.blast_prog == "blastn"
        assert stage.in_suffix == "_denoise.fasta"
        assert stage.out_suffix == "_blast.csv"
        assert stage.lineage_path is None

        expected_params = (
            "-db /path/to/db -max_target_seqs 1 -evalue 1e-05 "
            "-qcov_hsp_perc 90 -perc_identity 90 "
            '-outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sskingdoms" '
            "-num_threads 4"
        )
        assert stage.params == expected_params

    def test_init_custom(self, mock_config, tmp_dirs, lineage_file):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config,
            heading="custom_blast",
            blast_prog="blastn_custom",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_zotus.fa",
            out_suffix="_taxonomy.csv",
            db_path="/custom/db/path",
            lineage_path=lineage_file,
            maxhitnum=5,
            evalue=1e-10,
            qcov_hsp_perc=95,
            perc_identity=97,
            outfmt="6",
            specifiers="qseqid sseqid pident",
        )

        assert stage.heading == "custom_blast"
        assert stage.blast_prog == "blastn_custom"
        assert stage.in_suffix == "_zotus.fa"
        assert stage.out_suffix == "_taxonomy.csv"
        assert stage.lineage_path == lineage_file

        expected_params = (
            "-db /custom/db/path -max_target_seqs 5 -evalue 1e-10 "
            "-qcov_hsp_perc 95 -perc_identity 97 "
            '-outfmt "6 qseqid sseqid pident" -num_threads 4'
        )
        assert stage.params == expected_params

    def test_parse_params_remote_database(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir, db_path="nt"
        )

        assert "-remote" in stage.params
        assert "-num_threads" not in stage.params

    def test_parse_genus2otherlv(self, assign_taxa_stage, sample_lineage_data):
        expected_mapping = {
            "Salmo": [
                "Salmonidae",
                "Salmoniformes",
                "Actinopterygii",
                "Chordata",
                "Animalia",
            ],
            "Oncorhynchus": [
                "Salmonidae",
                "Salmoniformes",
                "Actinopterygii",
                "Chordata",
                "Animalia",
            ],
        }

        assert hasattr(assign_taxa_stage, "genus2otherlv")
        assert assign_taxa_stage.genus2otherlv == expected_mapping

    def test_parse_genus2otherlv_file_not_found(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        stage = AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            db_path="/path/to/db",
            lineage_path="/nonexistent/path.csv",
        )

        assert not hasattr(stage, "genus2otherlv")

    def test_setup(self, assign_taxa_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_denoise.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n")

        assign_taxa_stage.setup(test_prefix)

        assert len(assign_taxa_stage.runners) == 2
        blast_runner = assign_taxa_stage.runners[0]

        expected_infile = os.path.join(in_dir, f"{test_prefix}_denoise.fasta")
        expected_outfile = os.path.join(out_dir, f"{test_prefix}_blast.csv")

        expected_command = (
            f"blastn -query {expected_infile} "
            f"-db /path/to/database -max_target_seqs 1 -evalue 1e-05 "
            f"-qcov_hsp_perc 90 -perc_identity 90 "
            f'-outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sskingdoms" '
            f"-num_threads 4 -out {expected_outfile}"
        )

        assert blast_runner.command == expected_command
        assert blast_runner.prog_name == "Taxonomic assignment"

    def test_setup_no_lineage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AssignTaxaStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            db_path="/path/to/db",
        )

        input_file = os.path.join(in_dir, "test_denoise.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1\nACGT\n")

        stage.setup("test")

        assert len(stage.runners) == 1

    @patch("pandas.read_csv")
    def test_add_taxonomy(self, mock_read_csv, assign_taxa_stage, tmp_dirs):
        _in_dir, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assign_taxa_stage.blast_outfile = blast_file

        mock_blast_df = MagicMock()
        mock_blast_df.__getitem__.return_value = [
            "seq1|Salmo_trutta",
            "seq2|Oncorhynchus_mykiss",
        ]
        mock_read_csv.return_value = mock_blast_df

        with patch("numpy.array") as mock_array:
            mock_array.return_value.T = [
                ["Salmo_trutta", "Oncorhynchus_mykiss"],
                ["Salmo", "Oncorhynchus"],
                ["Salmonidae", "Salmonidae"],
                ["Salmoniformes", "Salmoniformes"],
                ["Actinopterygii", "Actinopterygii"],
                ["Chordata", "Chordata"],
                ["Animalia", "Animalia"],
            ]

            result = assign_taxa_stage.add_taxonomy()

            assert result is True
            mock_read_csv.assert_called_once_with(blast_file, header=None)
            mock_blast_df.to_csv.assert_called_once_with(
                blast_file, index=False, header=None
            )

    @patch("pandas.read_csv")
    def test_add_taxonomy_empty_result(
        self, mock_read_csv, assign_taxa_stage, tmp_dirs
    ):
        _, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assign_taxa_stage.blast_outfile = blast_file

        mock_read_csv.side_effect = pd.errors.EmptyDataError("No data")

        result = assign_taxa_stage.add_taxonomy()

        assert result is False
        assign_taxa_stage.config.logger.error.assert_called()

    def test_add_taxonomy_genus_not_found(self, assign_taxa_stage, tmp_dirs):
        _, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assign_taxa_stage.blast_outfile = blast_file

        mock_blast_data = MagicMock()
        mock_blast_data.__getitem__.return_value = ["seq1|Unknown_species"]

        with patch("pandas.read_csv", return_value=mock_blast_data):
            with patch("numpy.array") as mock_array:
                mock_array.return_value.T = []

                assign_taxa_stage.add_taxonomy()
                assign_taxa_stage.config.logger.error.assert_called()

    def test_add_taxonomy_family_level_identification(
        self, assign_taxa_stage, tmp_dirs
    ):
        _, out_dir = tmp_dirs
        blast_file = os.path.join(out_dir, "test_blast.csv")
        assign_taxa_stage.blast_outfile = blast_file

        assign_taxa_stage.genus2otherlv["TestGenus"] = [
            "Testidae",
            "TestOrder",
            "TestClass",
            "TestPhylum",
            "TestKingdom",
        ]

        mock_blast_data = MagicMock()
        mock_blast_data.__getitem__.return_value = ["seq1|Testidae_species"]

        with patch("pandas.read_csv", return_value=mock_blast_data):
            with patch("numpy.array") as mock_array:
                mock_array.return_value.T = []

                result = assign_taxa_stage.add_taxonomy()

                assert result is True
