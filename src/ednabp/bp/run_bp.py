import os
import re

from ..common import base_logger, config
from ..common.default_settings import SETTINGS
from .step_exec import (
    addhap,
    addlineage,
    blast,
    cutprimer,
    decompress,
    denoise,
    dereplicate,
    fqtofa,
    lca,
    merge,
)

STAGE_INPUT_SUFFIX = {
    "decompress": ".fastq.gz",
    "merge": "_R1.fastq",
    "cutprimer": ".fastq",
    "fqtofa": ".fastq",
    "dereplicate": ".fasta",
    "denoise": ".fasta",
    "blast": ".fasta",
    "addlineage": ".csv",
    "lca": ".csv",
    "addhap": ".csv",
}


class BioPipeline:
    def __init__(self, input_path: str, output_path: str, **settings):
        """
        A pipeline for processing eDNA bioinformatics workflows.

        :param input_path (str): The input path for raw sequences. This can be either a directory containing files or the path to a single file.
        :param output_path (str): The directory where output files will be saved.
        :param enabled_stages (list[str]): The process stages to be executed. The execution order will match the list order. Default:
         ["decompress", "merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "blast", "addlineage", "addhap"] (all stages will be run).
        :param settings: Additional optional arguments to configure pipeline stages and runtime behavior. These include:
          Input File Suffix:
            - raw_suffix (str): File suffix for raw input sequences. Default: "AUTO" (auto-detected based on starting stage).

          Merge Settings:
            - maxdiff (int): Maximum number of mismatches in the alignment. Default: 5.
            - pctid (int): Minimum %id of alignment. Default: 90.

          Cut Primer Settings:
            - rm_p_5 (str): Non-internal 5’ primer. Default: "GTCGGTAAAACTCGTGCCAGC" (MiFish-UF).
            - rm_p_3 (str): Non-internal 3’ primer. Default: "CAAACTGGGATTAGATACCCCACTATG" (reverse-complement MiFish-UR).
            - error_rate (float): The maximum rate of error could be tolerated. The actual error rate is computed as the number of errors in the match divided by the length of the matching part of the primer. Default: 0.15.
            - min_read_len (int): Discard processed reads that are shorter than this parameter. Default: 204.
            - max_read_len (int): Discard processed reads that are longer than this parameter. Default: 254.

          Denoise Settings:
            - minsize (int): Discard sequences with abundance that are smaller than this parameter. Default: 8.
            - alpha (int): Denoising sensitivity parameter. See UNOISE2 paper for definition. Default: 2.

          Blast Settings:
            - evalue (float): Expectation value (E) threshold for saving hits. Default: 0.00001.
            - qcov_hsp_perc (int): The %threshold of the query sequence that has to form an alignment against the reference to be retained. Default: 90.
            - perc_identity (int): Minimum percentage identity required for taxonomic assignment. Default: 90.
            - maxhitnum (int): Maximum number of hits to keep per query sequence. Default: 20.
            - specifiers (str): Output format specifiers for BLAST results. Default:
              "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".
            - blast_db (str): Path to the taxonomic database. Default: 'nt'.

          Add lineage Settings:
            - lineage_db (str): Path to the taxonomic lineage file. Default: 'nucleotide'.
            - entrez_email (str): The email used by NCBI to contact you in case of excessive usage or issues. Default: None

          LCA Settings:
            - tol_pct (float): Percentage of the top bitscore used as the inclusion threshold before LCA. Hits with bitscore >= top * (1 - tol_pct / 100) are included in the consensus. Default: 1.0.
            - score_column (str): Column used for score-based filtering. Default: 'bitscore'.
            - qseqid_column (str): Column used to group hits by query sequence. Default: 'qseqid'.

          External Program Setting:
            - usearch_prog (str): Command to execute USEARCH for merge, dereplicate, and denoise stages. Default: "usearch".
            - cutadapt_prog (str): Command to execute Cutadapt for cutprimer stage. Default: "cutadapt".
            - blast_prog (str): Command to execute BLAST for blast stage. Default: "blastn".

          Configuration Settings:
            - verbose (bool): Whether to enable verbose logging. Default: True.
            - dry (bool): If True, perform a dry run without executing commands. Default: False.
            - logger (Logger): Logger object for pipeline logging. Default: base_logger.logger.
            - n_cpu (int): Number of CPU cores to be used for processing. Default: 1.
            - memory (int): Maximum memory (in GB) allowed for processing. Default: 8.
        """
        assert os.path.exists(input_path), (
            f"Error: input path does not exist: {input_path}"
        )
        os.makedirs(output_path, exist_ok=True)
        input_is_dir = True if os.path.isdir(input_path) else False

        self.indir_path = (
            input_path if input_is_dir else os.path.dirname(input_path)
        )
        self.outdir_path = output_path

        self.add_default_settings(settings)

        self.add_config()

        self.determine_raw_suffix()

        self.setup_stages()

        if input_is_dir:
            self.run_stages_files()
        else:
            self.run_stages_one_file(input_path)

        self.close_file_handler()

    def add_default_settings(self, settings):
        self.stage_class = {
            "decompress": decompress.DecompressStage,
            "merge": merge.MergeStage,
            "cutprimer": cutprimer.CutPrimerStage,
            "fqtofa": fqtofa.FqToFaStage,
            "dereplicate": dereplicate.DereplicateStage,
            "denoise": denoise.DenoiseStage,
            "blast": blast.BlastStage,
            "addlineage": addlineage.AddLineageStage,
            "lca": lca.LcaStage,
            "addhap": addhap.AddHapStage,
        }

        self.enabled_stages = settings.get(
            "enabled_stages", SETTINGS["enabled_stages"]
        )
        self.stage_dir_name = SETTINGS["dir_name"].copy()
        self.stage_dir = {
            k: os.path.join(self.outdir_path, v)
            for k, v in self.stage_dir_name.items()
        }
        self.stage_suffix = SETTINGS["suffix"].copy()
        if "raw_suffix" in settings:
            self.stage_suffix["raw"] = settings["raw_suffix"]
        self.merge_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["merge"].items()
        }
        self.cutprimer_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["cutprimer"].items()
        }
        self.denoise_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["denoise"].items()
        }
        self.blast_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["blast"].items()
        }
        self.addlineage_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["addlineage"].items()
        }
        self.lca_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["lca"].items()
        }
        self.addhap_settings = {}
        self.prog_settings = {
            k: settings.get(f"{k}_prog", v)
            for k, v in SETTINGS["prog"].items()
        }
        self.config_basic_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["config_basic"].items()
        }
        self.config_machine_settings = {
            k: settings.get(k, v)
            for k, v in SETTINGS["config_machine"].items()
        }

    def add_config(self):
        self.config = config.Config(**self.config_basic_settings)
        self.config.add_machine_config(**self.config_machine_settings)
        fp_fh = base_logger.get_file_handler(
            os.path.join(self.outdir_path, "stages.log")
        )
        self.config.logger.addHandler(fp_fh)

    def determine_raw_suffix(self):
        raw_suffix = self.stage_suffix["raw"]
        if raw_suffix != SETTINGS["suffix"]["raw"]:
            self.config.logger.info(
                f"Using user-specified raw_suffix: {raw_suffix}"
            )
            print()
            return

        start_stage = self.enabled_stages[0]
        suffix = STAGE_INPUT_SUFFIX[start_stage]
        if "merge" in self.enabled_stages and start_stage == "decompress":
            suffix = f"_R1{suffix}"
        self.stage_suffix["raw"] = suffix
        self.config.logger.info(f"Using auto-determined raw_suffix: {suffix}")
        print()

    def setup_stages(self):
        self.stages = {}
        curr_dir = self.indir_path
        curr_suffix = self.stage_suffix["raw"]
        for stage in self.enabled_stages:
            if stage == "fqtofa":
                m = re.match(r".*\.(fastq|fq)$", curr_suffix)
                if m is not None:
                    self.stage_suffix["fqtofa"] = curr_suffix.replace(
                        m.group(1), "fasta"
                    )
                else:
                    self.config.logger.warning(
                        f"The in_suffix '{curr_suffix}' does not match the expected format (.fq/.fastq) for the 'fqtofa' stage."
                        "Skipping the stage"
                    )
                    continue

            stage_args = {
                "config": self.config,
                "in_dir": curr_dir,
                "out_dir": self.stage_dir[stage],
                "in_suffix": curr_suffix,
                "out_suffix": self.stage_suffix[stage],
            }

            if stage in ["merge", "dereplicate", "denoise"]:
                stage_args["usearch_prog"] = self.prog_settings["usearch"]
            elif stage == "cutprimer":
                stage_args["cutadapt_prog"] = self.prog_settings["cutadapt"]
            elif stage in ["blast", "assigntaxa"]:
                stage_args["blast_prog"] = self.prog_settings["blast"]

            if stage == "addhap":
                stage_args["denoise_dir"] = self.stage_dir["denoise"]

            if stage in [
                "merge",
                "cutprimer",
                "denoise",
                "blast",
                "addlineage",
                "lca",
            ]:
                stage_args.update(eval(f"self.{stage}_settings"))

            self.stages[stage] = self.stage_class[stage](**stage_args)

            curr_dir = self.stage_dir[stage]
            curr_suffix = self.stage_suffix[stage]

    def run_one_file(self, prefix):
        base_logger.logger.info(f"Sample ID: {prefix}")
        for k, s in self.stages.items():
            s.setup(prefix)
            is_complete = s.run()
            if not is_complete:
                self.config.logger.error(
                    f"Error: process errors at stage: {k}\n"
                )
                break
            if self.config.verbose:
                print()

    def run_stages_files(self):
        suffix = self.stage_suffix["raw"]
        files = os.listdir(self.indir_path)
        prefixes = [
            file.replace(suffix, "") for file in files if file.endswith(suffix)
        ]

        if len(prefixes) == 0:
            self.config.logger.warning(
                f"No input files found with suffix '{suffix}' in directory '{self.indir_path}'"
            )
            return

        for prefix in prefixes:
            self.run_one_file(prefix)

    def run_stages_one_file(self, input_path):
        prefix = os.path.basename(input_path).replace(
            self.stage_suffix["raw"], ""
        )
        self.run_one_file(prefix)

    def close_file_handler(self):
        base_logger.close_file_handler(self.config.logger)
