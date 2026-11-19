import os
import shutil

import pytest

from ednabp.bp.step_exec.qfilter import QFilterStage, qscore_to_error_rate
from ednabp.bp.step_exec.utils import read_fasta

USEARCH = shutil.which("usearch")


@pytest.mark.parametrize(
    "qscore, rate", [(10, 0.1), (20, 0.01), (30, 0.001), (0, 1.0)]
)
def test_qscore_to_error_rate(qscore, rate):
    assert qscore_to_error_rate(qscore) == pytest.approx(rate)


class TestQFilterStage:
    @pytest.fixture
    def stage(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        return QFilterStage(config=mock_config, in_dir=in_dir, out_dir=out_dir)

    def test_init(self, stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        assert stage.in_dir == in_dir
        assert stage.out_dir == out_dir
        assert stage.usearch_prog == "usearch"
        assert stage.in_suffix == ".fastq"
        assert stage.out_suffix == ".fasta"
        # defaults from the settings: Q30, scores up to 70
        assert stage.maxee_rate == pytest.approx(0.001)
        assert stage.params == "-fastq_maxee_rate 0.001 -fastq_qmax 70"

    def test_init_custom(self, mock_config, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        stage = QFilterStage(
            config=mock_config,
            usearch_prog="usearch12",
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix=".fq",
            min_qscore=15,
            max_qscore=50,
            maxns=0,
        )
        assert stage.usearch_prog == "usearch12"
        assert stage.in_suffix == ".fq"
        assert stage.params == (
            "-fastq_maxee_rate 0.0316228 -fastq_qmax 50 -fastq_maxns 0"
        )

    def test_setup(self, stage, tmp_dirs):
        in_dir, out_dir = tmp_dirs
        infile = os.path.join(in_dir, "s1.fastq")
        with open(infile, "w") as f:
            f.write("@r1\nACGT\n+\nIIII\n")

        stage.setup("s1")

        outfile = os.path.join(out_dir, "s1.fasta")
        assert len(stage.runners) == 2
        assert stage.runners[0].prog_name == "Filter low-quality reads"
        assert stage.runners[0].command == (
            f"usearch -fastq_filter {infile} -fastq_maxee_rate 0.001"
            f" -fastq_qmax 70 -fastaout {outfile}"
        )
        assert stage.runners[1].runner is stage.runners[0]

    def test_setup_missing_infile(self, stage):
        with pytest.raises(FileNotFoundError):
            stage.setup("does_not_exist")


@pytest.mark.skipif(USEARCH is None, reason="usearch must be on PATH")
def test_run_filters_low_quality(mock_config, tmp_dirs):
    """min_qscore 23 is an error rate of 0.005: 'lowq' (Q20, 0.01 per base)
    is dropped, 'hasN' by maxns 0; 'ontQ50' ('S') is above usearch's
    default limit of 42 but within max_qscore 50."""
    in_dir, out_dir = tmp_dirs
    with open(os.path.join(in_dir, "s1.fastq"), "w") as f:
        f.write(
            "@good\nACGTACGTAC\n+\nIIIIIIIIII\n"
            "@lowq\nACGTACGTAC\n+\n5555555555\n"
            "@hasN\nACGTNCGTAC\n+\nIIIIIIIIII\n"
            "@ontQ50\nACGTACGTAC\n+\nSSSSSSSSSS\n"
        )
    stage = QFilterStage(
        config=mock_config,
        usearch_prog=USEARCH,
        in_dir=in_dir,
        out_dir=out_dir,
        min_qscore=23,
        max_qscore=50,
        maxns=0,
    )
    stage.setup("s1")
    assert stage.run()

    records = list(read_fasta(os.path.join(out_dir, "s1.fasta")))
    assert records == [("good", "ACGTACGTAC"), ("ontQ50", "ACGTACGTAC")]


@pytest.mark.skipif(USEARCH is None, reason="usearch must be on PATH")
def test_score_above_max_qscore_fails(mock_config, tmp_dirs):
    in_dir, out_dir = tmp_dirs
    with open(os.path.join(in_dir, "s1.fastq"), "w") as f:
        f.write("@ontQ50\nACGTACGTAC\n+\nSSSSSSSSSS\n")
    stage = QFilterStage(
        config=mock_config,
        usearch_prog=USEARCH,
        in_dir=in_dir,
        out_dir=out_dir,
        max_qscore=49,
    )
    stage.setup("s1")
    assert not stage.run()
