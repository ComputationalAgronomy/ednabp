import os

import pandas as pd
import pytest

from ednabp.bp.step_exec.lca import (
    RANKS,
    LcaStage,
    filter_by_score,
    lca_of_lineages,
    normalize_species_names,
)


def blast_row(
    qseqid,
    bitscore,
    kingdom="Eukaryota",
    phylum="Chordata",
    cls="Actinopteri",
    order="Gadiformes",
    family="Gadidae",
    genus="Gadus",
    species="Gadus_morhua",
    pident=99.0,
):
    return {
        "qseqid": qseqid,
        "bitscore": bitscore,
        "pident": pident,
        "kingdom": kingdom,
        "phylum": phylum,
        "class": cls,
        "order": order,
        "family": family,
        "genus": genus,
        "species": species,
    }


class TestLcaStage:
    @pytest.fixture
    def lca_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return LcaStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

    def test_normalize_simple_binomial(self):
        assert normalize_species_names("Homo_sapiens") == ["Homo_sapiens"]

    def test_normalize_subspecies_returns_full_then_base(self):
        assert normalize_species_names("Homo_sapiens_neanderthalensis") == [
            "Homo_sapiens_neanderthalensis",
            "Homo_sapiens",
        ]

    def test_normalize_hybrid_splits_on_x(self):
        assert normalize_species_names("Quercus_x_rosacea") == [
            "Quercus",
            "rosacea",
        ]

    def test_normalize_sp_notation(self):
        assert "Gadus_sp." in normalize_species_names("Gadus_sp.-A12345")

    def test_normalize_empty_string(self):
        assert normalize_species_names("") == []

    def test_filter_tol_zero_keeps_only_best(self):
        g = pd.DataFrame({"bitscore": [300, 280, 250], "qseqid": ["q1"] * 3})
        assert list(
            filter_by_score(g, "bitscore", tol_pct=0.0)["bitscore"]
        ) == [300]

    def test_filter_tol_nonzero_keeps_within_window(self):
        g = pd.DataFrame({"bitscore": [300, 295, 280], "qseqid": ["q1"] * 3})
        assert set(
            filter_by_score(g, "bitscore", tol_pct=2.0)["bitscore"]
        ) == {
            300,
            295,
        }

    def test_filter_empty_group(self):
        g = pd.DataFrame({"bitscore": [], "qseqid": []})
        assert filter_by_score(g, "bitscore", tol_pct=0.0).empty

    def test_filter_all_equal_scores_keeps_all(self):
        g = pd.DataFrame({"bitscore": [200, 200, 200], "qseqid": ["q1"] * 3})
        assert len(filter_by_score(g, "bitscore", tol_pct=0.0)) == 3

    def test_lca_identical_lineages(self):
        lineage = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        lca, rank = lca_of_lineages([lineage, lineage])
        assert rank == "species"
        assert lca[RANKS.index("species")] == "Gadus_morhua"

    def test_lca_diverge_at_species(self):
        l1 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        l2 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_chalcogrammus",
        ]
        lca, rank = lca_of_lineages([l1, l2])
        assert rank == "genus"
        assert lca[RANKS.index("species")] == ""

    def test_lca_diverge_at_genus(self):
        l1 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        l2 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Melanogrammus",
            "Melanogrammus_aeglefinus",
        ]
        lca, rank = lca_of_lineages([l1, l2])
        assert rank == "family"
        assert lca[RANKS.index("genus")] == ""

    def test_lca_diverge_at_kingdom(self):
        l1 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        l2 = ["Bacteria", "Proteobacteria", "", "", "", "", ""]
        lca, rank = lca_of_lineages([l1, l2])
        assert rank == "NA"
        assert all(v == "" for v in lca)

    def test_lca_empty(self):
        lca, rank = lca_of_lineages([])
        assert rank == "NA"
        assert lca == [""] * len(RANKS)

    def test_lca_subspecies_with_species(self):
        l1 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua_callarias",
        ]
        l2 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        lca, rank = lca_of_lineages([l1, l2])
        assert rank == "species"
        assert lca[RANKS.index("species")] == "Gadus_morhua"

    def test_lca_single_lineage(self):
        lineage = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        lca, rank = lca_of_lineages([lineage])
        assert rank == "species"
        assert lca == [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]

    def test_lca_missing_middle_rank(self):
        l1 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        l2 = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "",
            "Gadidae",
            "Gadus",
            "Gadus_morhua",
        ]
        lca, rank = lca_of_lineages([l1, l2])
        assert rank == "species"
        assert lca[RANKS.index("order")] == ""
        assert lca[RANKS.index("family")] == "Gadidae"
        assert lca[RANKS.index("genus")] == "Gadus"
        assert lca[RANKS.index("species")] == "Gadus_morhua"

    def test_lca_intraspecific_assignment(self):
        lineage = [
            "Eukaryota",
            "Chordata",
            "Actinopteri",
            "Gadiformes",
            "Gadidae",
            "Gadus",
            "Gadus_morhua_callarias",
        ]
        lca, rank = lca_of_lineages([lineage, lineage])
        assert rank == "species"
        assert lca[RANKS.index("species")] == "Gadus_morhua_callarias"

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = LcaStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.in_suffix == ".csv"
        assert stage.out_suffix == ".csv"
        assert stage.tol_pct == 1.0
        assert stage.score_column == "bitscore"
        assert stage.qseqid_column == "qseqid"
        assert os.path.isdir(out_dir)

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = LcaStage(
            config=mock_config,
            heading="custom_lca",
            in_dir=in_dir,
            out_dir=out_dir,
            tol_pct=5.0,
            score_column="evalue",
            qseqid_column="query_id",
        )

        assert stage.heading == "custom_lca"
        assert stage.tol_pct == 5.0
        assert stage.score_column == "evalue"
        assert stage.qseqid_column == "query_id"

    def test_setup(self, lca_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"
        input_file = os.path.join(in_dir, f"{test_prefix}.csv")
        pd.DataFrame([blast_row("Zotu1", 300)]).to_csv(input_file, index=False)

        lca_stage.setup(test_prefix)

        assert len(lca_stage.runners) == 1
        assert lca_stage.infile == input_file
        assert lca_stage.outfile == os.path.join(out_dir, f"{test_prefix}.csv")

    def test_setup_missing_infile_raises(self, lca_stage):
        with pytest.raises(FileNotFoundError):
            lca_stage.setup("nonexistent")

    def test_setup_dry_skips_runner(self, dry_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = LcaStage(config=dry_config, in_dir=in_dir, out_dir=out_dir)
        input_file = os.path.join(in_dir, "sample001.csv")
        pd.DataFrame([blast_row("Zotu1", 300)]).to_csv(input_file, index=False)

        stage.setup("sample001")

        assert len(stage.runners) == 0

    def test_run_lca_single_hit_resolves_to_species(self, lca_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        pd.DataFrame([blast_row("Zotu1", 300)]).to_csv(
            lca_stage.infile, index=False
        )

        result = lca_stage.run_lca()

        assert result is True
        output = pd.read_csv(lca_stage.outfile)
        assert (
            output.loc[output["qseqid"] == "Zotu1", "lca_rank"].iloc[0]
            == "species"
        )

    def test_run_lca_diverging_species_resolves_to_genus(
        self, lca_stage, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        rows = [
            blast_row("Zotu1", 300, species="Gadus_morhua"),
            blast_row("Zotu1", 300, species="Gadus_chalcogrammus"),
        ]
        pd.DataFrame(rows).to_csv(lca_stage.infile, index=False)

        lca_stage.run_lca()

        output = pd.read_csv(lca_stage.outfile)
        row = output.loc[output["qseqid"] == "Zotu1"].iloc[0]
        assert row["lca_rank"] == "genus"
        assert pd.isna(row["species"])

    def test_run_lca_tight_tol_pct_excludes_low_score_hit(
        self, lca_stage, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        rows = [
            blast_row("Zotu1", 300, species="Gadus_morhua"),
            blast_row("Zotu1", 300, species="Gadus_morhua"),
            blast_row(
                "Zotu1",
                250,
                genus="Melanogrammus",
                species="Melanogrammus_aeglefinus",
            ),
        ]
        pd.DataFrame(rows).to_csv(lca_stage.infile, index=False)

        lca_stage.run_lca()

        output = pd.read_csv(lca_stage.outfile)
        assert (
            output.loc[output["qseqid"] == "Zotu1", "lca_rank"].iloc[0]
            == "species"
        )

    def test_run_lca_wide_tol_pct_causes_family_lca(
        self, mock_config, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        stage = LcaStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir, tol_pct=2.0
        )
        stage.infile = os.path.join(in_dir, "test.csv")
        stage.outfile = os.path.join(out_dir, "test.csv")
        rows = [
            blast_row("Zotu1", 300, species="Gadus_morhua"),
            blast_row(
                "Zotu1",
                295,
                genus="Melanogrammus",
                species="Melanogrammus_aeglefinus",
            ),
        ]
        pd.DataFrame(rows).to_csv(stage.infile, index=False)

        stage.run_lca()

        output = pd.read_csv(stage.outfile)
        assert (
            output.loc[output["qseqid"] == "Zotu1", "lca_rank"].iloc[0]
            == "family"
        )

    def test_run_lca_multiple_queries(self, lca_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        rows = [
            blast_row("Zotu1", 300, species="Gadus_morhua"),
            blast_row("Zotu2", 300, species="Gadus_chalcogrammus"),
        ]
        pd.DataFrame(rows).to_csv(lca_stage.infile, index=False)

        lca_stage.run_lca()

        output = pd.read_csv(lca_stage.outfile)
        assert len(output) == 2
        assert set(output["qseqid"]) == {"Zotu1", "Zotu2"}

    def test_run_lca_num_assignments_counts_filtered_hits(
        self, lca_stage, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        rows = [
            blast_row("Zotu1", 300),
            blast_row("Zotu1", 300),
            blast_row("Zotu1", 100),
        ]
        pd.DataFrame(rows).to_csv(lca_stage.infile, index=False)

        lca_stage.run_lca()

        output = pd.read_csv(lca_stage.outfile)
        assert (
            output.loc[output["qseqid"] == "Zotu1", "num_assignments"].iloc[0]
            == 2
        )

    def test_run_lca_best_row_metadata_from_highest_score(
        self, mock_config, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        stage = LcaStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir, tol_pct=7.0
        )
        stage.infile = os.path.join(in_dir, "test.csv")
        stage.outfile = os.path.join(out_dir, "test.csv")
        rows = [
            blast_row("Zotu1", 300, pident=99.0),
            blast_row("Zotu1", 280, pident=85.0),
        ]
        pd.DataFrame(rows).to_csv(stage.infile, index=False)

        stage.run_lca()

        output = pd.read_csv(stage.outfile)
        assert (
            output.loc[output["qseqid"] == "Zotu1", "pident"].iloc[0] == 99.0
        )

    def test_run_lca_output_has_expected_columns(self, lca_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        pd.DataFrame([blast_row("Zotu1", 300)]).to_csv(
            lca_stage.infile, index=False
        )

        lca_stage.run_lca()

        output = pd.read_csv(lca_stage.outfile)
        assert "lca_rank" in output.columns
        assert "num_assignments" in output.columns

    def test_run_lca_returns_false_on_missing_taxonomy_columns(
        self, lca_stage, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs
        lca_stage.infile = os.path.join(in_dir, "test.csv")
        lca_stage.outfile = os.path.join(out_dir, "test.csv")
        pd.DataFrame({"qseqid": ["Zotu1"], "bitscore": [300]}).to_csv(
            lca_stage.infile, index=False
        )

        result = lca_stage.run_lca()

        assert result is False
        lca_stage.config.logger.error.assert_called()
