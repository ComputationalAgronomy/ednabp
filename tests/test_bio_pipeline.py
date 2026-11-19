import os
import tempfile
from unittest.mock import Mock, patch

import pytest

from ednabp.bp.run_bp import BioPipeline
from ednabp.common.config import Config
from ednabp.common.default_settings import SETTINGS


def mock_add_config(self):
    self.config = Mock()
    self.config.logger = Mock()
    self.config.logger.handlers = []


class TestBioPipeline:
    @pytest.fixture
    def sample_fastq_files(self, tmp_dirs):
        in_dir, _ = tmp_dirs

        sample_files = [
            "sample001_R1.fastq.gz",
            "sample001_R2.fastq.gz",
            "sample002_R1.fastq.gz",
            "sample002_R2.fastq.gz",
        ]

        for filename in sample_files:
            filepath = os.path.join(in_dir, filename)
            with open(filepath, "w") as f:
                f.write("@read1\nACGT\n+\nIIII\n")

        return in_dir

    def test_init_input_dir(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
        ):
            pipeline = BioPipeline(in_dir, out_dir, verbose=True, dry=True)

            assert pipeline.indir_path == in_dir
            assert pipeline.outdir_path == out_dir
            assert os.path.exists(out_dir)

    def test_init_input_single_file(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "single_sample_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("@read1\nACGT\n+\nIIII\n")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_one_file"),
        ):
            pipeline = BioPipeline(input_file, out_dir, verbose=True, dry=True)

            assert pipeline.indir_path == in_dir
            assert pipeline.outdir_path == out_dir

    def test_init_input_nonexistent_path(self, tmp_dirs):
        _in_dir, out_dir = tmp_dirs
        nonexistent_path = "/nonexistent/path"

        with pytest.raises(
            FileNotFoundError, match="input path does not exist"
        ):
            BioPipeline(nonexistent_path, out_dir)

    def test_add_default_settings(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        custom_settings = {
            "enabled_stages": ["decompress", "merge"],
            "maxdiff": 10,
            "verbose": False,
            "n_cpu": 8,
        }

        with (
            patch.object(BioPipeline, "add_config"),
            patch.object(BioPipeline, "determine_raw_suffix"),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
        ):
            pipeline = BioPipeline(in_dir, out_dir, **custom_settings)
            pipeline.config = Mock()

            assert pipeline.enabled_stages == ["decompress", "merge"]
            assert pipeline.stage_dir_name["decompress"] == "decompress"
            assert pipeline.merge_settings["maxdiff"] == 10
            assert pipeline.config_basic_settings["verbose"] is False
            assert pipeline.config_machine_settings["n_cpu"] == 8

    def test_add_config(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch("ednabp.common.config.Config") as mock_config_class,
            patch(
                "ednabp.common.base_logger.get_file_handler"
            ) as mock_file_handler,
            patch.object(BioPipeline, "close_file_handler"),
        ):
            mock_config = Mock(spec=Config)
            mock_config.logger = Mock()
            mock_config_class.return_value = mock_config
            mock_handler = Mock()
            mock_file_handler.return_value = mock_handler

            _pipeline = BioPipeline(in_dir, out_dir, verbose=True)

            mock_config_class.assert_called_once()
            mock_config.add_machine_config.assert_called_once()
            mock_config.logger.addHandler.assert_called_once_with(mock_handler)

            expected_log_path = os.path.join(out_dir, "stages.log")
            mock_file_handler.assert_called_once_with(expected_log_path)

    def test_determine_raw_suffix_user_specified(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test.custom_suffix")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["fqtofa"],
                raw_suffix=".custom_suffix",
                verbose=False,
            )

            assert pipeline.stage_suffix["raw"] == ".custom_suffix"

    def test_determine_raw_suffix_decompress_no_merge(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch("ednabp.common.base_logger.get_file_handler"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["decompress"],
                verbose=False,
            )

            assert pipeline.stage_suffix["raw"] == ".fastq.gz"

    def test_determine_raw_suffix_decompress_with_merge(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch("ednabp.common.base_logger.get_file_handler"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["decompress", "merge", "cutprimer"],
                verbose=False,
            )

            assert pipeline.stage_suffix["raw"] == "_R1.fastq.gz"

    def test_determine_raw_suffix_merge_start(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch("ednabp.common.base_logger.get_file_handler"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["merge", "cutprimer"],
                verbose=False,
            )

            assert pipeline.stage_suffix["raw"] == "_R1.fastq"

    def test_determine_raw_suffix_fqtofa_no_merge(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test.fastq")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch("ednabp.common.base_logger.get_file_handler"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["fqtofa", "dereplicate"],
                verbose=False,
            )

            assert pipeline.stage_suffix["raw"] == ".fastq"

    def test_determine_raw_suffix_dereplicate_start(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test.fasta")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch("ednabp.common.base_logger.get_file_handler"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["dereplicate", "denoise", "blast"],
                verbose=False,
            )

            assert pipeline.stage_suffix["raw"] == ".fasta"

    @patch("ednabp.bp.step_exec.decompress.DecompressStage")
    @patch("ednabp.bp.step_exec.merge.MergeStage")
    def test_setup_stages(
        self, mock_merge_stage, mock_decompress_stage, tmp_dirs
    ):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["decompress", "merge"],
                verbose=True,
                dry=True,
            )

            mock_decompress_stage.assert_called_once()
            mock_merge_stage.assert_called_once()

            assert "decompress" in pipeline.stages
            assert "merge" in pipeline.stages

    def test_setup_stages_fqtofa_suffix_handling(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test.fastq")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.fqtofa.FqToFaStage"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["fqtofa"],
            )

            expected_fqtofa_suffix = ".fasta"
            assert pipeline.stage_suffix["fqtofa"] == expected_fqtofa_suffix

    def test_default_stages_use_qfilter(self):
        stages = SETTINGS["enabled_stages"]
        assert "qfilter" in stages and "fqtofa" not in stages
        assert stages.index("cutprimer") < stages.index("qfilter")
        assert stages.index("qfilter") < stages.index("dereplicate")

    def test_setup_stages_qfilter(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        with open(os.path.join(in_dir, "test.fq"), "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "add_config", mock_add_config),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["qfilter", "dereplicate"],
                raw_suffix=".fq",
                min_qscore=30,
                max_qscore=50,
                usearch_prog="usearch12",
            )

        stage = pipeline.stages["qfilter"]
        assert pipeline.stage_suffix["qfilter"] == ".fasta"
        assert stage.in_suffix == ".fq" and stage.out_suffix == ".fasta"
        assert stage.params == "-fastq_maxee_rate 0.001 -fastq_qmax 50"
        assert stage.usearch_prog == "usearch12"
        assert pipeline.stages["dereplicate"].in_suffix == ".fasta"

    def test_setup_stages_fqtofa_invalid_suffix(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fasta")
        with open(input_file, "w") as f:
            f.write("test")

        def mock_determine_suffix(self):
            self.stage_suffix["raw"] = ".fasta"

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch.object(
                BioPipeline, "determine_raw_suffix", mock_determine_suffix
            ),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["fqtofa"],
            )

            assert "fqtofa" not in pipeline.stages

    @patch("ednabp.common.base_logger.logger")
    def test_run_one_file(self, mock_logger, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "add_config"),
            patch.object(BioPipeline, "determine_raw_suffix"),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
        ):
            pipeline = BioPipeline(in_dir, out_dir, verbose=True, dry=True)
            pipeline.config = Mock()

            mock_stage1 = Mock()
            mock_stage1.run.return_value = True
            mock_stage2 = Mock()
            mock_stage2.run.return_value = True

            pipeline.stages = {"stage1": mock_stage1, "stage2": mock_stage2}

            pipeline.run_stages(["test_prefix"])

            mock_stage1.setup.assert_called_once_with("test_prefix")
            mock_stage1.run.assert_called_once()
            mock_stage2.setup.assert_called_once_with("test_prefix")
            mock_stage2.run.assert_called_once()

            mock_logger.info.assert_called_with("Sample ID: test_prefix")

    def test_run_one_file_stage_failure(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "add_config"),
            patch.object(BioPipeline, "determine_raw_suffix"),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
        ):
            pipeline = BioPipeline(in_dir, out_dir, verbose=True, dry=True)
            pipeline.config = Mock()
            pipeline.config.logger = Mock()

            mock_stage1 = Mock()
            mock_stage1.run.return_value = True
            mock_stage2 = Mock()
            mock_stage2.run.return_value = False
            mock_stage3 = Mock()

            pipeline.stages = {
                "stage1": mock_stage1,
                "stage2": mock_stage2,
                "stage3": mock_stage3,
            }

            pipeline.run_stages(["test_prefix"])

            mock_stage1.run.assert_called_once()
            mock_stage2.run.assert_called_once()
            mock_stage3.run.assert_not_called()

            pipeline.config.logger.error.assert_called()

    def test_run_stages_files(self, sample_fastq_files, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with (
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages") as mock_run_stages,
            patch.object(BioPipeline, "close_file_handler"),
        ):
            _pipeline = BioPipeline(
                in_dir,
                out_dir,
                verbose=True,
                dry=True,
            )

            # All samples go to the stages together, in sorted order
            mock_run_stages.assert_called_once_with(["sample001", "sample002"])

    def test_run_stages_one_file(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "single_sample_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_stages") as mock_run_stages,
            patch.object(BioPipeline, "close_file_handler"),
        ):
            _pipeline = BioPipeline(
                input_file,
                out_dir,
                verbose=True,
                dry=True,
            )

            mock_run_stages.assert_called_once_with(["single_sample"])

    def test_stage_parameter_passing(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        custom_settings = {
            "enabled_stages": ["merge", "cutprimer", "denoise", "blast"],
            "maxdiff": 15,
            "pctid": 85,
            "rm_p_5": "CUSTOM5PRIMER",
            "rm_p_3": "CUSTOM3PRIMER",
            "minsize": 10,
            "alpha": 3,
            "evalue": 1e-10,
            "blast_db": "/custom/db/path",
            "task": "dc-megablast",
            "usearch_prog": "custom_usearch",
            "cutadapt_prog": "custom_cutadapt",
            "blast_prog": "custom_blastn",
        }

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.merge.MergeStage") as mock_merge,
            patch(
                "ednabp.bp.step_exec.cutprimer.CutPrimerStage"
            ) as mock_cutprimer,
            patch("ednabp.bp.step_exec.denoise.DenoiseStage") as mock_denoise,
            patch("ednabp.bp.step_exec.blast.BlastStage") as mock_assigntaxa,
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                **custom_settings,
                verbose=True,
                dry=True,
            )
            pipeline.config = Mock()

            merge_call_kwargs = mock_merge.call_args[1]
            assert merge_call_kwargs["maxdiff"] == 15
            assert merge_call_kwargs["pctid"] == 85
            assert merge_call_kwargs["usearch_prog"] == "custom_usearch"

            cutprimer_call_kwargs = mock_cutprimer.call_args[1]
            assert cutprimer_call_kwargs["rm_p_5"] == "CUSTOM5PRIMER"
            assert cutprimer_call_kwargs["rm_p_3"] == "CUSTOM3PRIMER"
            assert cutprimer_call_kwargs["cutadapt_prog"] == "custom_cutadapt"

            denoise_call_kwargs = mock_denoise.call_args[1]
            assert denoise_call_kwargs["minsize"] == 10
            assert denoise_call_kwargs["alpha"] == 3
            assert denoise_call_kwargs["usearch_prog"] == "custom_usearch"

            assigntaxa_call_kwargs = mock_assigntaxa.call_args[1]
            assert assigntaxa_call_kwargs["evalue"] == 1e-10
            assert assigntaxa_call_kwargs["blast_db"] == "/custom/db/path"
            assert assigntaxa_call_kwargs["task"] == "dc-megablast"
            assert assigntaxa_call_kwargs["blast_prog"] == "custom_blastn"

    def test_cluster_stage_parameter_passing(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with open(os.path.join(in_dir, "test.fasta"), "w") as f:
            f.write(">seq1;size=10\nACGT\n")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.denoise.DenoiseStage") as mock_denoise,
            patch("ednabp.bp.step_exec.cluster.ClusterStage") as mock_cluster,
        ):
            BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["cluster"],
                minsize=10,
                cluster_minsize=4,
                usearch_prog="custom_usearch",
            )

            mock_denoise.assert_not_called()
            cluster_call_kwargs = mock_cluster.call_args[1]
            assert cluster_call_kwargs["minsize"] == 4
            assert "cluster_minsize" not in cluster_call_kwargs
            assert cluster_call_kwargs["usearch_prog"] == "custom_usearch"
            assert cluster_call_kwargs["in_suffix"] == ".fasta"
            assert cluster_call_kwargs["out_dir"] == os.path.join(
                out_dir, "cluster"
            )

    def test_otutab_parameter_passing(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with open(os.path.join(in_dir, "test.fasta"), "w") as f:
            f.write(">seq1;size=10\nACGT\n")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.denoise.DenoiseStage") as mock_denoise,
            patch("ednabp.bp.step_exec.cluster.ClusterStage") as mock_cluster,
            patch(
                "ednabp.bp.step_exec.cluster_ont.ClusterOntStage"
            ) as mock_ont,
        ):
            BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["denoise", "cluster", "cluster_ont"],
                otutab_id=0.99,
                otutab_strand=False,
            )

            for mock_stage in (mock_denoise, mock_cluster, mock_ont):
                kwargs = mock_stage.call_args[1]
                assert kwargs["otutab_id"] == 0.99
                assert kwargs["otutab_strand"] is False

    def test_otutab_id_unset_uses_stage_defaults(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with open(os.path.join(in_dir, "test.fasta"), "w") as f:
            f.write(">seq1;size=10\nACGT\n")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
        ):
            pipeline = BioPipeline(
                in_dir, out_dir, enabled_stages=["denoise", "cluster_ont"]
            )

        assert pipeline.stages["denoise"].otutab_id == 0.97
        assert pipeline.stages["cluster_ont"].otutab_id == 0.90

    def test_cluster_ont_parameter_passing(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with open(os.path.join(in_dir, "test.fasta"), "w") as f:
            f.write(">read1\nACGT\n")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch(
                "ednabp.bp.step_exec.cluster_ont.ClusterOntStage"
            ) as mock_stage,
        ):
            BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["cluster_ont"],
                cluster_id=0.88,
                medaka=True,
                medaka_model="r1041_e82_400bps_sup_v5.0.0",
                cluster_strand=True,
                strand=False,
                racon_prog="/opt/racon",
                usearch_prog="/opt/usearch",
            )

            kwargs = mock_stage.call_args[1]
            assert kwargs["cluster_id"] == 0.88
            assert kwargs["medaka"] is True
            assert kwargs["medaka_model"] == "r1041_e82_400bps_sup_v5.0.0"
            assert kwargs["cluster_strand"] is True
            assert "strand" not in kwargs
            assert "otutab_id" not in kwargs
            assert kwargs["otutab_strand"] is True
            assert kwargs["racon_prog"] == "/opt/racon"
            assert kwargs["usearch_prog"] == "/opt/usearch"
            assert kwargs["out_dir"] == os.path.join(out_dir, "cluster_ont")

    @pytest.mark.parametrize(
        "stages, expected",
        [
            (["denoise", "blast", "addhap"], "denoise"),
            (["cluster", "blast", "addhap"], "cluster"),
            (["cluster_ont", "blast", "addhap"], "cluster_ont"),
            (["blast", "addhap"], None),
        ],
    )
    def test_addhap_reads_enabled_otu_stage(self, tmp_dirs, stages, expected):
        """Without an OTU stage, a run starting at blast reads the OTU
        files from its input directory (expected=None)."""
        in_dir, out_dir = tmp_dirs

        with open(os.path.join(in_dir, "test.fasta"), "w") as f:
            f.write(">seq1;size=10\nACGT\n")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.addhap.AddHapStage") as mock_addhap,
        ):
            BioPipeline(in_dir, out_dir, enabled_stages=stages)

            expected_dir = (
                in_dir if expected is None else os.path.join(out_dir, expected)
            )
            assert mock_addhap.call_args[1]["otu_dir"] == expected_dir

    def test_stage_directory_creation(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.decompress.DecompressStage"),
            patch("ednabp.bp.step_exec.merge.MergeStage"),
        ):
            pipeline = BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["decompress", "merge"],
                verbose=True,
                dry=True,
            )

            expected_decompress_dir = os.path.join(out_dir, "decompress")
            expected_merge_dir = os.path.join(out_dir, "merge")

            assert pipeline.stage_dir["decompress"] == expected_decompress_dir
            assert pipeline.stage_dir["merge"] == expected_merge_dir

    def test_eval_usage_in_stage_args(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch("ednabp.bp.step_exec.merge.MergeStage") as mock_merge,
        ):
            BioPipeline(
                in_dir,
                out_dir,
                enabled_stages=["merge"],
                verbose=True,
                dry=True,
            )

            mock_merge.assert_called_once()

    def test_out_directory_creation(self, tmp_dirs):
        in_dir, _ = tmp_dirs

        input_file = os.path.join(in_dir, "test_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with tempfile.TemporaryDirectory() as temp_base:
            new_out_dir = os.path.join(temp_base, "new_out_dir")

            assert not os.path.exists(new_out_dir)

            with (
                patch.object(BioPipeline, "add_config"),
                patch.object(BioPipeline, "determine_raw_suffix"),
                patch.object(BioPipeline, "setup_stages"),
                patch.object(BioPipeline, "run_stages_files"),
                patch.object(BioPipeline, "close_file_handler"),
            ):
                BioPipeline(in_dir, new_out_dir, verbose=True, dry=True)

                assert os.path.exists(new_out_dir)


class FakeStage:
    """Records setup calls; fails for the (stage, sample) pairs in `fail`."""

    def __init__(self, name, calls, fail=()):
        self.name, self.calls, self.fail = name, calls, fail
        self.runners = []

    def setup(self, prefix):
        self.prefix = prefix
        self.calls.append((self.name, prefix))

    def run(self):
        return (self.name, self.prefix) not in self.fail


class FakePool(FakeStage):
    pools_samples = True

    def run_all(self, prefixes):
        self.calls.append((self.name, tuple(prefixes)))
        return ["_pooled"]


class TestPooling:
    def make_pipeline(self, tmp_dirs, **settings):
        in_dir, out_dir = tmp_dirs
        for name in ["B.fasta", "A.fasta", "A_otutab_query.fasta"]:
            with open(os.path.join(in_dir, name), "w") as f:
                f.write(">Uniq1;size=10;\nACGT\n")
        with (
            patch.object(BioPipeline, "run_stages_files"),
            patch.object(BioPipeline, "close_file_handler"),
            patch.object(BioPipeline, "add_config", mock_add_config),
        ):
            return BioPipeline(in_dir, out_dir, **settings)

    def fake_stages(self, pipeline, fail=()):
        calls = []
        pipeline.stages = {
            k: (FakePool if k == "pool" else FakeStage)(k, calls, fail)
            for k in pipeline.stages
        }
        return calls

    def test_stage_by_stage(self, tmp_dirs):
        """Each stage runs for all samples before the next; OTU side files
        in the input folder are not samples."""
        pipeline = self.make_pipeline(
            tmp_dirs, enabled_stages=["dereplicate", "denoise", "blast"]
        )
        calls = self.fake_stages(pipeline)

        pipeline.run_stages_files()

        assert calls == [
            ("dereplicate", "A"),
            ("dereplicate", "B"),
            ("denoise", "A"),
            ("denoise", "B"),
            ("blast", "A"),
            ("blast", "B"),
        ]

    @pytest.mark.parametrize(
        "stages",
        [
            ["dereplicate", "pool", "denoise", "blast", "addhap"],
            ["dereplicate", "denoise", "pool", "blast", "addhap"],
        ],
    )
    def test_pool_before_or_after_denoise(self, tmp_dirs, stages):
        pipeline = self.make_pipeline(tmp_dirs, enabled_stages=stages)
        calls = self.fake_stages(pipeline)

        pipeline.run_stages_files()

        at = stages.index("pool")
        before = [(s, p) for s in stages[:at] for p in ("A", "B")]
        after = [(s, "_pooled") for s in stages[at + 1 :]]
        assert calls == before + [("pool", ("A", "B"))] + after

    def test_failed_sample_not_pooled(self, tmp_dirs):
        pipeline = self.make_pipeline(
            tmp_dirs,
            enabled_stages=["dereplicate", "denoise", "pool", "blast"],
        )
        calls = self.fake_stages(pipeline, fail={("denoise", "B")})

        pipeline.run_stages_files()

        assert ("pool", ("A",)) in calls
        assert ("blast", "_pooled") in calls
        pipeline.config.logger.error.assert_called()

    def test_setup_error_drops_sample(self, tmp_dirs):
        pipeline = self.make_pipeline(
            tmp_dirs, enabled_stages=["dereplicate", "denoise"]
        )
        calls = self.fake_stages(pipeline)

        def missing(prefix):
            calls.append(("dereplicate", prefix))
            if prefix == "A":
                raise FileNotFoundError("A.fasta not found")

        pipeline.stages["dereplicate"].setup = missing
        pipeline.stages["dereplicate"].prefix = "B"

        pipeline.run_stages(["A", "B"])

        assert ("denoise", "A") not in calls
        assert ("denoise", "B") in calls

    @pytest.mark.parametrize(
        "stages, expected",
        [
            (["dereplicate", "pool", "denoise", "blast", "addhap"], "denoise"),
            (["dereplicate", "denoise", "pool", "blast", "addhap"], "pool"),
        ],
    )
    def test_addhap_reads_blast_input(self, tmp_dirs, stages, expected):
        pipeline = self.make_pipeline(tmp_dirs, enabled_stages=stages)
        _, out_dir = tmp_dirs
        assert pipeline.stages["addhap"].otu_dir == os.path.join(
            out_dir, expected
        )

    def test_pool_name_setting(self, tmp_dirs):
        pipeline = self.make_pipeline(
            tmp_dirs,
            enabled_stages=["dereplicate", "pool", "denoise"],
            pool_name="farm",
        )
        assert pipeline.stages["pool"].pool_name == "farm"
        _, out_dir = tmp_dirs
        assert pipeline.stages["denoise"].in_dir == os.path.join(
            out_dir, "pool"
        )


def test_dry_run_sets_up_every_stage(tmp_dirs):
    """A dry run writes no files, yet every stage is set up for every
    sample and each command is logged instead of run."""
    from unittest.mock import Mock

    from ednabp.common.config import Config

    in_dir, out_dir = tmp_dirs
    for s in ["A", "B"]:
        for r in ["R1", "R2"]:
            with open(os.path.join(in_dir, f"{s}_{r}.fastq"), "w") as f:
                f.write("@r\nACGT\n+\nIIII\n")
    logger = Mock()

    def add_config(self):
        self.config = Config(logger=logger, verbose=True, dry=True)
        self.config.add_machine_config(n_cpu=1, memory=1)

    with (
        patch.object(BioPipeline, "add_config", add_config),
        patch.object(BioPipeline, "close_file_handler"),
    ):
        BioPipeline(
            in_dir,
            out_dir,
            enabled_stages=[
                "merge",
                "cutprimer",
                "qfilter",
                "dereplicate",
                "pool",
                "denoise",
                "blast",
                "addlineage",
                "lca",
                "addhap",
            ],
            lineage_db=None,
            dry=True,
        )
    logged = [c.args[0] for c in logger.info.call_args_list]
    dry_runs = [m for m in logged if m.startswith("DRY RUN: ")]
    assert not logger.error.called, logger.error.call_args_list
    # one merge per sample, then the pool and the stages on the pool
    assert sum("-fastq_mergepairs" in m for m in dry_runs) == 2
    assert any("-unoise3" in m and "_pooled" in m for m in dry_runs)
    assert any(m.startswith("DRY RUN: blastn") for m in dry_runs)
    written = [
        f
        for _, _, files in os.walk(out_dir)
        for f in files
        if f != "stages.log"
    ]
    assert written == []
