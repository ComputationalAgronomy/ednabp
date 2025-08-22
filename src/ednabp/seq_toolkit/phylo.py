import os

from ..bp.step_build import stage_builder
from ..common import base_logger, config
from ..data import BPData, MitoData


class Phylo:
    def __init__(
        self,
        out_dir,
        iqtree_prog="iqtree2",
        verbose: bool = True,
        n_cpu: int = "AUTO",
        model: str = "TEST",
        bootstrap: int = None,
        overwrite: bool = False,
    ):
        self.out_dir = out_dir
        self.IQTREE_PROG = iqtree_prog
        self.add_config(verbose, n_cpu, model, bootstrap, overwrite)

    def add_config(self, verbose, n_cpu, model, bootstrap, overwrite):
        self.config = config.Config()
        self.config.verbose = verbose
        self.config.add_machine_info(n_cpu=n_cpu)
        self.config.add_iqtree_config(model, bootstrap, overwrite)
        fp_fh = base_logger.get_file_handler(
            os.path.join(self.out_dir, "phylo.log")
        )
        self.config.logger.addHandler(fp_fh)

    def reconstruct(self, data: str | BPData | MitoData):
        if isinstance(data, BPData | MitoData):
            from tempfile import TemporaryDirectory

            from ..write import seq_writer

            writer = seq_writer.SeqWriter(
                data, self.config.n_cpu, self.config.verbose
            )
            with TemporaryDirectory() as tmp_dir:
                tmp_fasta = os.path.join(tmp_dir, "tmp.fasta")
                derep_fasta = os.path.join(tmp_dir, "mltree.fasta")
                in_fasta = os.path.join(self.out_dir, "mltree.aln")
                writer.fasta(tmp_fasta)
                writer.derep_fasta(tmp_fasta, derep_fasta)
                writer.align_fasta(derep_fasta, in_fasta)
        elif isinstance(data, str):
            if not os.path.exists(data):
                raise FileNotFoundError(f"data not found: {data}")
            in_fasta = data
        else:
            raise ValueError(
                "data should be a BPData, MitoData class, or valid .FASTA file."
            )

        in_dir = os.path.dirname(in_fasta)
        in_basename = os.path.basename(in_fasta)
        prefix, in_suffix = os.path.splitext(in_basename)
        stage = PhyloReconStage(
            self.config,
            iqtree_prog=self.IQTREE_PROG,
            in_dir=in_dir,
            out_dir=self.out_dir,
            in_suffix=in_suffix,
        )

        stage.setup(prefix)
        stage.run()


class PhyloReconStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        iqtree_prog="iqtree2",
        in_dir="",
        out_dir="",
        in_suffix=".fasta",
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.IQTREE_PROG = iqtree_prog
        self.in_suffix = in_suffix
        self.parse_params()

    def parse_params(
        self,
    ):
        config = self.config.get_iqtree_config()
        self.params = f"-m {config['model']} -nt {config['threads']}"
        if config["bootstrap"]:
            self.params += f" -b {config['bootstrap']}"
        if config["overwrite"]:
            self.params += " -redo"

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        phylo_outfile = os.path.join(self.out_dir, prefix)
        cmd = (
            f"{self.IQTREE_PROG} -s {self.infile}"
            f" --prefix {phylo_outfile} {self.params}"
        )
        super().add_stage("Reconstruct phylogenetic tree", cmd)

    def run(self):
        super().run()
        return all(self.output)
