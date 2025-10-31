import csv
import os
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from ednabp.bp.step_exec.addlineage import AddLineageStage, get_empty_rank


class TestAddLineageStage:
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
    def addlineage_stage(self, mock_config, tmp_dirs, lineage_file):
        in_dir, out_dir = tmp_dirs
        return AddLineageStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            lineage_db=lineage_file,
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AddLineageStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.in_suffix == "_blast.csv"
        assert stage.out_suffix == "_taxa.csv"
        assert stage.lineage_db == "nucleotide"

    def test_init_custom(self, mock_config, tmp_dirs, lineage_file):
        in_dir, out_dir = tmp_dirs
        stage = AddLineageStage(
            config=mock_config,
            heading="custom_lineage",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_blast_results.csv",
            out_suffix="_taxonomy.csv",
            lineage_db=lineage_file,
        )

        assert stage.heading == "custom_lineage"
        assert stage.in_suffix == "_blast_results.csv"
        assert stage.out_suffix == "_taxonomy.csv"
        assert stage.lineage_db == lineage_file

    def test_parse_genus2otherlv(self, addlineage_stage, sample_lineage_data):
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

        assert hasattr(addlineage_stage, "genus2otherlv")
        assert addlineage_stage.genus2otherlv == expected_mapping

    def test_parse_lineage_db_invalid_file(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with pytest.raises(ValueError):
            AddLineageStage(
                config=mock_config,
                in_dir=in_dir,
                out_dir=out_dir,
                lineage_db="/nonexistent/path.csv",
            )

    def test_setup(self, addlineage_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}_blast.csv")
        with open(input_file, "w") as f:
            f.write("sseqid,pident\nseq1|Salmo_trutta,95\n")

        addlineage_stage.setup(test_prefix)

        assert len(addlineage_stage.runners) >= 1

    def test_setup_no_lineage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        addlineage_stage = AddLineageStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            lineage_db=None,
        )

        input_file = os.path.join(in_dir, "test_blast.csv")
        with open(input_file, "w") as f:
            f.write("sseqid,pident\nseq1,95\n")

        addlineage_stage.setup("test")

        assert len(addlineage_stage.runners) == 0

    @patch("pandas.read_csv")
    def test_get_lineage_from_custom_db(
        self, mock_read_csv, addlineage_stage, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        blast_file = os.path.join(in_dir, "test_blast.csv")
        output_file = os.path.join(out_dir, "test_taxa.csv")
        addlineage_stage.infile = blast_file
        addlineage_stage.outfile = output_file

        mock_blast_df = pd.DataFrame(
            {"sseqid": ["seq1|Salmo_trutta", "seq2|Unknown_species"]}
        )
        mock_read_csv.return_value = mock_blast_df

        with patch.object(mock_blast_df, "merge") as mock_merge:
            mock_merged = MagicMock()
            mock_merge.return_value = mock_merged

            result = addlineage_stage.get_lineage_from_custom_db()

            assert result is True
            mock_merged.to_csv.assert_called_once_with(output_file, index=False)

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
        stage = AddLineageStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            lineage_db="nucleotide",
        )

        assert stage.lineage_db == "nucleotide"
        assert hasattr(stage, "get_lineage_func")
        assert stage.entrez_email == "your.email@example.com"