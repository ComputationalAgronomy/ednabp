import os

import pandas as pd
import pytest

from ednabp.bp.step_exec.addhap import AddHapStage


def write_inputs(in_dir, denoise_dir, prefix, taxa, fasta, otutab):
    with open(os.path.join(in_dir, f"{prefix}.csv"), "w") as f:
        f.write(taxa)
    with open(os.path.join(denoise_dir, f"{prefix}.fasta"), "w") as f:
        f.write(fasta)
    with open(os.path.join(denoise_dir, f"{prefix}_otu_table.txt"), "w") as f:
        f.write(otutab)


class TestAddHapStage:
    @pytest.fixture
    def denoise_dir(self, tmp_dirs):
        _, out_dir = tmp_dirs
        path = os.path.join(out_dir, "denoise")
        os.makedirs(path, exist_ok=True)
        return path

    @pytest.fixture
    def addhap_stage(self, mock_config, tmp_dirs, denoise_dir):
        in_dir, out_dir = tmp_dirs
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
        assert stage.otutab_suffix == "_otu_table.txt"

    def test_setup(self, addhap_stage, tmp_dirs, denoise_dir):
        in_dir, _ = tmp_dirs
        write_inputs(
            in_dir,
            denoise_dir,
            "sample001",
            "qseqid,sseqid\nZotu1,seq1\n",
            ">Zotu1\nACGT\n>Zotu2\nACGG\n",
            "#OTU ID\tsample001\nZotu1\t87\nZotu2\t42\n",
        )

        addhap_stage.setup("sample001")

        assert len(addhap_stage.runners) == 1
        assert addhap_stage.otutab_file == os.path.join(
            denoise_dir, "sample001_otu_table.txt"
        )

    def test_setup_missing_otu_table(
        self, addhap_stage, tmp_dirs, denoise_dir
    ):
        in_dir, _ = tmp_dirs
        write_inputs(
            in_dir, denoise_dir, "s1", "qseqid\nZotu1\n", ">Zotu1\nA\n", ""
        )
        os.remove(os.path.join(denoise_dir, "s1_otu_table.txt"))

        with pytest.raises(FileNotFoundError, match="s1_otu_table.txt"):
            addhap_stage.setup("s1")

    def test_add_haplotype_info(self, addhap_stage, tmp_dirs, denoise_dir):
        in_dir, out_dir = tmp_dirs
        write_inputs(
            in_dir,
            denoise_dir,
            "test",
            "qseqid,sseqid\nZotu1,gb|EF173208.1|\n",
            ">Zotu1\nACGT\n>Zotu2\nACGG\n",
            "#OTU ID\ttest\nZotu2\t42\nZotu1\t87\n",
        )
        addhap_stage.setup("test")

        assert addhap_stage.add_haplotype_info() is True

        output_df = pd.read_csv(os.path.join(out_dir, "test.csv"))
        assert output_df["qseqid"].tolist() == ["Zotu1", "Zotu2"]
        assert output_df["zotu"].tolist() == ["ACGT", "ACGG"]
        assert output_df["size"].tolist() == [87, 42]
        assert output_df["size"].dtype == "int64"

    def test_add_haplotype_info_cluster_otus(
        self, addhap_stage, tmp_dirs, denoise_dir
    ):
        """OTU labels from the cluster/cluster_ont stages ("Otu<n>")."""
        in_dir, out_dir = tmp_dirs
        write_inputs(
            in_dir,
            denoise_dir,
            "s1",
            "qseqid,sseqid\nOtu2,hit\n",
            ">Otu1\nAAAA\n>Otu2\nCCCC\n>Otu10\nGGGG\n",
            "#OTU ID\ts1\nOtu1\t120\nOtu2\t59\nOtu10\t3\n",
        )
        addhap_stage.setup("s1")
        addhap_stage.add_haplotype_info()

        output_df = pd.read_csv(os.path.join(out_dir, "s1.csv"))
        # Row with a BLAST hit first, then the others by OTU number
        assert output_df["qseqid"].tolist() == ["Otu2", "Otu1", "Otu10"]
        assert output_df["size"].tolist() == [59, 120, 3]

    def test_read_otu_table(self, addhap_stage, tmp_dirs):
        _, out_dir = tmp_dirs
        path = os.path.join(out_dir, "t_otu_table.txt")
        with open(path, "w") as f:
            f.write("#OTU ID\ts1\ts2\nZotu1\t87\t3\nZotu2\t0\t2\n")

        assert addhap_stage.read_otu_table(path) == {"Zotu1": 90, "Zotu2": 2}

    def test_read_otu_table_bad_header(self, addhap_stage, tmp_dirs):
        _, out_dir = tmp_dirs
        path = os.path.join(out_dir, "t.txt")
        with open(path, "w") as f:
            f.write("Uniq1;size=87;\tchfilter zotu\n")

        with pytest.raises(ValueError):
            addhap_stage.read_otu_table(path)

    def test_read_fasta_sequences(self, addhap_stage, tmp_dirs):
        _, out_dir = tmp_dirs
        fasta_file = os.path.join(out_dir, "test.fasta")

        with open(fasta_file, "w") as f:
            f.write(">Zotu1\nACGT\n>Zotu2\nTGCA\n")

        sequences = addhap_stage.read_fasta_sequences(fasta_file)

        assert sequences == {"Zotu1": "ACGT", "Zotu2": "TGCA"}
