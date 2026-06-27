import os
import tempfile
from unittest.mock import Mock, patch

import pytest

from ednabp.bp.run_bp import BioPipeline
from ednabp.common.config import Config


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

        with pytest.raises(AssertionError, match="input path does not exist"):
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

            pipeline.run_one_file("test_prefix")

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

            pipeline.run_one_file("test_prefix")

            mock_stage1.run.assert_called_once()
            mock_stage2.run.assert_called_once()
            mock_stage3.run.assert_not_called()

            pipeline.config.logger.error.assert_called()

    def test_run_stages_files(self, sample_fastq_files, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        with (
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_one_file") as mock_run_one_file,
            patch.object(BioPipeline, "close_file_handler"),
        ):
            _pipeline = BioPipeline(
                in_dir,
                out_dir,
                verbose=True,
                dry=True,
            )

            expected_samples = {"sample001", "sample002"}

            actual_calls = mock_run_one_file.call_args_list[:2]
            actual_samples = {call[0][0] for call in actual_calls}
            assert actual_samples == expected_samples

    def test_run_stages_one_file(self, tmp_dirs):
        in_dir, out_dir = tmp_dirs

        input_file = os.path.join(in_dir, "single_sample_R1.fastq.gz")
        with open(input_file, "w") as f:
            f.write("test")

        with (
            patch.object(BioPipeline, "add_config", mock_add_config),
            patch.object(BioPipeline, "setup_stages"),
            patch.object(BioPipeline, "run_one_file") as mock_run_one_file,
            patch.object(BioPipeline, "close_file_handler"),
        ):
            _pipeline = BioPipeline(
                input_file,
                out_dir,
                verbose=True,
                dry=True,
            )

            mock_run_one_file.assert_called_once()
            call_args = mock_run_one_file.call_args[0]
            assert call_args[0] == "single_sample"

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
            assert assigntaxa_call_kwargs["blast_prog"] == "custom_blastn"

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
            pipeline = BioPipeline(
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
