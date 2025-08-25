import os
from abc import ABC
from collections import defaultdict

from ..bp.step_build import stage_builder
from ..bp.step_exec import dereplicate
from ..common import base_writer


def filter_by_dict(target_dict, filter_dict):
    match = True
    for filter_key, filter_value in filter_dict.items():
        if filter_key not in filter_dict:
            continue

        if isinstance(filter_value, list):
            if target_dict[filter_key] not in filter_value:
                match = False
                break
        else:
            if target_dict[filter_key] != filter_value:
                match = False
                break

    return match


class Writer(base_writer.BaseWriter, ABC):
    """
    An abstract class for running sequence related analysis
    """

    def __init__(self, data, verbose=False, n_cpu=1):
        super().__init__(data, verbose)
        self.config.add_machine_config(n_cpu=n_cpu)

    def load_units2fasta(
        self, add_taxa_info: bool = True, filter_taxa=None, filter_sample=None
    ):
        self.units2fasta = defaultdict(str)
        for sample_id in self.sample_id_used:
            for hap, lv_dict in self.data.sample_data[
                sample_id
            ].hap2level.items():
                if isinstance(filter_taxa, dict):
                    if not filter_by_dict(lv_dict, filter_taxa):
                        continue
                if isinstance(filter_sample, dict):
                    if not filter_by_dict(
                        self.data.sample_metadata[sample_id], filter_sample
                    ):
                        continue
                unit_name = (
                    f"{lv_dict['class']}_{lv_dict['family']}_{lv_dict['species']}-"
                    if add_taxa_info
                    else ""
                )
                title = f"{unit_name}{sample_id}_{hap}"
                seq = self.data.sample_data[sample_id].hap_seq[hap]
                self.units2fasta[unit_name] += f">{title}\n{seq}\n"

    def fasta(self, out_path: str, add_taxa_info: bool = True):
        """
        Write sequences to a FASTA file.

        :param save_path: Path to the output FASTA file.
        """
        self.load_units2fasta(add_taxa_info)
        fasta_str = "".join(list(self.units2fasta.values()))
        with open(out_path, "w") as file:
            file.write(fasta_str)
        self.config.logger.info(f"Saved sequence file to: {out_path}")

    def derep_fasta(self, in_path, out_path):
        in_dir = os.path.dirname(in_path)
        out_dir = os.path.dirname(out_path)
        in_basename = os.path.basename(in_path)
        out_basename = os.path.basename(out_path)
        stage = dereplicate.DereplicateStage(
            self.config,
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix=in_basename,
            out_suffix=out_basename,
            annot_size=False,
            write_report=False,
        )
        stage.setup("")
        stage.run()

    def align_fasta(self, in_path, out_path):
        in_dir = os.path.dirname(in_path)
        out_dir = os.path.dirname(out_path)
        in_basename = os.path.basename(in_path)
        out_basename = os.path.basename(out_path)
        stage = MSAStage(
            self.config,
            in_dir=in_dir,
            out_dir=out_dir,
            in_suffix=in_basename,
            out_suffix=out_basename,
        )
        stage.setup("")
        stage.run()


class MSAStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        clustalo_prog="clustalo",
        in_dir="",
        out_dir="",
        in_suffix=".fasta",
        out_suffix="_aligned.fasta",
        overwrite=True,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.CLUSTALO_PROG = clustalo_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.parse_params(overwrite)

    def parse_params(self, overwrite):
        self.params = "--force" if overwrite else ""

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        msa_outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.check_infile()
        self.check_outdir()
        cmd = (
            f"{self.CLUSTALO_PROG} -i {self.infile}"
            f" -o {msa_outfile} {self.params}"
        )
        super().add_stage("Align sequences", cmd)

    def run(self):
        super().run()
        return all(self.output)
