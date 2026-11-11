import os

import pytest

from ednabp.bp.step_exec.cluster import ClusterStage


class TestClusterStage:
    @pytest.fixture
    def cluster_stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return ClusterStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

    def test_init(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = ClusterStage(
            config=mock_config, in_dir=in_dir, out_dir=out_dir
        )

        assert stage.config == mock_config
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.usearch_prog == "usearch"
        assert stage.in_suffix == ".fasta"
        assert stage.out_suffix == ".fasta"
        assert stage.cluster_report_suffix == "_cluster_report.txt"
        assert stage.report_suffix == "_report.txt"
        assert stage.params == "-minsize 2 -relabel Otu"

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = ClusterStage(
            config=mock_config,
            heading="custom_cluster",
            usearch_prog="usearch12",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix="_unique.fa",
            out_suffix="_otus.fa",
            minsize=4,
            seq_label="",
        )

        assert stage.heading == "custom_cluster"
        assert stage.usearch_prog == "usearch12"
        assert stage.in_suffix == "_unique.fa"
        assert stage.out_suffix == "_otus.fa"
        assert stage.params == "-minsize 4"

    def test_setup(self, cluster_stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        test_prefix = "sample001"

        input_file = os.path.join(in_dir, f"{test_prefix}.fasta")
        with open(input_file, "w") as f:
            f.write(">seq1;size=10\nACGT\n>seq2;size=5\nTGCA\n")

        cluster_stage.setup(test_prefix)

        # Command + output redirection, otutab + report, row summing
        assert len(cluster_stage.runners) == 5
        command_runner = cluster_stage.runners[0]

        expected_outfile = os.path.join(out_dir, f"{test_prefix}.fasta")
        expected_report = os.path.join(
            out_dir, f"{test_prefix}_cluster_report.txt"
        )
        expected_command = (
            f"usearch -cluster_otus {input_file} "
            f"-minsize 2 -relabel Otu "
            f"-otus {expected_outfile} -uparseout {expected_report}"
        )

        assert command_runner.command == expected_command
        assert command_runner.prog_name == "Cluster unique sequences into OTUs"

        otutab_runner = cluster_stage.runners[2]
        expected_otutab_command = (
            f"usearch -otutab "
            f"{input_file} "
            f"-otus {expected_outfile} -id 0.97 -strand both -threads 4 "
            f"-otutabout {os.path.join(out_dir, f'{test_prefix}_otu_table.txt')} "
            f"-mapout {os.path.join(out_dir, f'{test_prefix}_otu_mapping.txt')} "
            f"-notmatched {os.path.join(out_dir, f'{test_prefix}_otu_notmatched.fasta')}"
        )
        assert otutab_runner.command == expected_otutab_command

    def test_setup_missing_infile(self, cluster_stage):
        with pytest.raises(FileNotFoundError):
            cluster_stage.setup("does_not_exist")
