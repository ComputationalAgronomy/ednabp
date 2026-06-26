import os
from unittest.mock import patch

import pandas as pd
import pytest

from ednabp.bp.step_exec.addhap import AddHapStage


class TestAddHapStage:
    @pytest.fixture
    def addhap_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        denoise_dir = os.path.join(out_dir, "denoise")
        os.makedirs(denoise_dir, exist_ok=True)
        return AddHapStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            denoise_dir=denoise_dir,
        )

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = AddHapStage(
            config=mock_config,
            in_dir=in_dir,
            out_dir=out_dir,
            denoise_dir="denoise",
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.in_suffix == ".csv"
        assert stage.out_suffix == ".csv"
        assert stage.denoise_dir == "denoise"

    def test_setup(self, addhap_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"
        taxa_file = os.path.join(in_dir, f"{test_prefix}.csv")
        addhap_stage.denoise_file = os.path.join(
            out_dir, "denoise", f"{test_prefix}.fasta"
        )
        addhap_stage.report_file = os.path.join(
            out_dir, "denoise", f"{test_prefix}_denoise_report.txt"
        )

        with open(taxa_file, "w") as f:
            f.write("qseqid,sseqid\nZotu1,seq1\n")

        with open(addhap_stage.denoise_file, "w") as f:
            f.write(">Zotu1\nACGT\n>Zotu2\nACGG\n")

        with open(addhap_stage.report_file, "w") as f:
            f.write("Uniq1;size=87;\tchfilter zotu\n")
            f.write("Uniq2;size=42;\tchfilter zotu\n")

        addhap_stage.setup(test_prefix)

        assert len(addhap_stage.runners) == 1

    def test_add_haplotype_info(self, addhap_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        taxa_file = os.path.join(in_dir, "test.csv")
        hap_file = os.path.join(out_dir, "test.csv")
        addhap_stage.infile = taxa_file
        addhap_stage.outfile = hap_file
        addhap_stage.denoise_file = os.path.join(
            out_dir, "denoise", "test.fasta"
        )
        addhap_stage.report_file = os.path.join(
            out_dir, "denoise", "test_denoise_report.txt"
        )

        with open(taxa_file, "w") as f:
            f.write("qseqid,sseqid\nZotu1,gb|EF173208.1|\n")

        with open(addhap_stage.denoise_file, "w") as f:
            f.write(">Zotu1\nACGT\n>Zotu2\nACGG\n")

        with open(addhap_stage.report_file, "w") as f:
            f.write("Uniq1;size=87;\tchfilter zotu\n")
            f.write("Uniq2;size=42;\tchfilter zotu\n")

        result = addhap_stage.add_haplotype_info()

        assert result is True
        assert os.path.exists(hap_file)

        output_df = pd.read_csv(hap_file)
        assert len(output_df) == 2
        assert "Zotu1" in output_df["qseqid"].values
        assert "Zotu2" in output_df["qseqid"].values

    def test_read_denoise_report(self, addhap_stage, tmp_dirs):
        _, out_dir = tmp_dirs
        report_file = os.path.join(out_dir, "test_denoise_report.txt")

        with open(report_file, "w") as f:
            f.write("Uniq1;size=87;\tchfilter zotu\n")
            f.write("Uniq2;size=42;\tchfilter chimera\n")
            f.write("Uniq3;size=2;\tchfilter zotu\n")

        zotu_sizes = addhap_stage.read_denoise_report(report_file)

        assert zotu_sizes == {"Zotu1": 87, "Zotu2": 2}

    def test_read_fasta_sequences(self, addhap_stage, tmp_dirs):
        _, out_dir = tmp_dirs
        fasta_file = os.path.join(out_dir, "test.fasta")

        with open(fasta_file, "w") as f:
            f.write(">Zotu1\nACGT\n>Zotu2\nTGCA\n")

        sequences = addhap_stage.read_fasta_sequences(fasta_file)

        assert sequences == {"Zotu1": "ACGT", "Zotu2": "TGCA"}
