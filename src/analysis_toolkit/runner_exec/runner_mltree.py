import os
import subprocess
import tempfile

from analysis_toolkit.runner_build import (base_runner, utils, utils_sequence)


class MLTreeRunner(base_runner.SequenceRunner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)

    @base_runner.log_execution("Reconstruct ML tree", "write_mltree.log")
    def run_write(self,
            target_list: list[str],
            target_level: str,
            unit_level:str = "species",
            save_dir:str = '.',
            save_prefix: str = "ml_tree",
            model: str = None,
            bootstrap: int = None,
            threads: int = None,
            dereplicate_sequence: bool = True,
            n_unit_threshold:int = 1,
            sample_id_list: list[str] = []
        ) -> None:
        """
        Reconstruct a phylogenetic tree for a list of targets using IQTREE and write a .TREEFILE file.
        (IQTREE2 command reference: http://www.iqtree.org/doc/Command-Reference)

        :param target_list: A list of targets to be plotted (e.g., ["FamilyA", "FamilyB", etc]).
        :param target_level: The taxonomic level of the targets (e.g., family, genus, species).
        :param units_level: The taxonomic level of the units. Default is "species"
        :param save_dir: Directory to save the output files.
        :param save_prefix: Prefix for the output file names. Default is "ml_tree".
        :param model: Model to specify for tree inference. If not specified, it will use the best-fit model found. Default is None.
        :param bootstrap: Number of bootstrap replicates. Default is None.
        :param threads: Number of threads to use. If not specified, it will automatically determine the best number of cores given the current data and computer. Default is None.
        :param dereplicate_sequence: Whether to dereplicate the sequence. Default is True.
        :param n_unit_threshold: Minimum number of seqeunces for an unit to be included in the analysis. Default is 1.
        :param sample_id_list: A list of sample IDs to plot. Default is None (plot all samples).
        """
        os.makedirs(save_dir, exist_ok=True)

        ml_fasta_path = os.path.join(save_dir, f'{save_prefix}.aln')

        self._load_sample_id_list(sample_id_list)

        self._write_mltree_fasta(
            target_list=target_list,
            target_level=target_level,
            unit_level=unit_level,
            save_path=ml_fasta_path,
            n_unit_threshold=n_unit_threshold,
            dereplicate_sequence=dereplicate_sequence,
        )

        self._run_iqtree2(
            seq_path=ml_fasta_path,
            save_dir=save_dir,
            save_prefix=save_prefix,
            model=model,
            bootstrap=bootstrap,
            threads=threads
        )

        self.analysis_type = "mltree_run"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "target_list": target_list,
                "target_level": target_level,
                "unit_level": unit_level,
                "save_prefix": save_prefix,
                "model": model,
                "bootstrap": bootstrap,
                "threads": threads,
                "dereplicate_sequence": dereplicate_sequence,
                "n_unit_threshold": n_unit_threshold,
            }
        )

    def run_plot(self):
        return super().run_plot()

    def _check_mltree_overwrite(self, save_dir: str, save_prefix: str) -> str:
        # TODO(SW): Eventually, replace with argparse.
        """
        Check if an MLTree run with the given prefix already exists in the specified directory.
        If it does, prompt the user to decide whether to redo the run or not.

        :param save_dir: Directory where the output files are saved.
        :param prefix: Prefix for the output file names.
        :return: User input choice ('-redo', '--redo-tree', '--undo', or 'stop').
        """
        ckp_path = os.path.join(save_dir, save_prefix + '.ckp.gz') # e.g. save/dir/SpA.ckp.gz
        if os.path.exists(ckp_path):
            print(
                f"> MLTree checkpoint fileCheckpoint ({ckp_path}) indicates that a previous run successfully finished  already exists.\n"
                "Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n"
                "Use `--redo-tree` option if you want to restore ModelFinder and only redo tree search.\n"
                "Use `--undo` option if you want to continue previous run when changing/adding options.\n"
            )
            user_input_choice = ['-redo', '--redo-tree', '--undo', 'stop']
            help_msg = f"({'/'.join(user_input_choice)}): "
            while True:
                user_input = input(help_msg).strip().lower()
                if user_input in user_input_choice:
                    return user_input
                print("> Invalid input.")
        else:
            pass

    def _write_mltree_fasta(self,
            target_list: list[str],
            target_level: str,
            unit_level: str,
            save_path: str,
            n_unit_threshold: int,
            dereplicate_sequence: bool,
        ) -> None:
        """
        Write an aligned FASTA file for a list of targets.
        """
        self.logger.info(f"Writing MLTree FASTA file: {save_path}...")
 
        try:
            temp_dir = tempfile.TemporaryDirectory()
            fasta_path = os.path.join(temp_dir.name, 'mltree.fa')

            for target_name in target_list:
                self._load_units2fasta_dict(
                    target_name=target_name,
                    target_level=target_level,
                    unit_level=unit_level,
                    n_unit_threshold=n_unit_threshold
                )

            utils_sequence.write_fasta(self.units2fasta, save_path=fasta_path, dereplicate=dereplicate_sequence)
            utils_sequence.align_fasta(seq_file=fasta_path, aln_file=save_path)

        finally:
            temp_dir.cleanup()

    def _run_iqtree2(self,
            seq_path:str,
            save_dir: str,
            save_prefix: str,
            model:str = None,
            bootstrap: int = None,
            threads: int = None
        ) -> None:
        """
        Run IQTREE2 command.
        """
        checkpoint = self._check_mltree_overwrite(save_dir, save_prefix)
        if checkpoint == 'stop':
            self.logger.info("Stopping the run.")
            return

        model = model or 'TEST'
        prefix_path = os.path.join(save_dir, save_prefix)
        threads = threads or 'AUTO'

        cmd = [
            'iqtree2', '-m', model, '-s', seq_path, '--prefix', prefix_path, '-nt', threads
        ]
        if bootstrap:
            cmd.extend(["-b", str(bootstrap)])
        if checkpoint:
            cmd.append(checkpoint)

        utils.run_subprocess("IQTREE2", cmd, save_dir)