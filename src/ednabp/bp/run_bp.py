import os
import re

from ..common import base_logger, config
from ..common.default_settings import SETTINGS
from .step_exec import (
    addhap,
    addlineage,
    blast,
    cluster,
    cluster_ont,
    cutprimer,
    decompress,
    denoise,
    dereplicate,
    fqtofa,
    lca,
    merge,
    qfilter,
)

STAGE_INPUT_SUFFIX = {
    "decompress": ".fastq.gz",
    "merge": "_R1.fastq",
    "cutprimer": ".fastq",
    "fqtofa": ".fastq",
    "qfilter": ".fastq",
    "dereplicate": ".fasta",
    "denoise": ".fasta",
    "cluster": ".fasta",
    "cluster_ont": ".fasta",
    "blast": ".fasta",
    "addlineage": ".csv",
    "lca": ".csv",
    "addhap": ".csv",
}


SIDE_FILE_SUFFIXES = [
    SETTINGS["suffix"]["members"],
    SETTINGS["suffix"]["otutab_query"],
]


class BioPipeline:
    def __init__(self, input_path: str, output_path: str, **settings):
        """
        A pipeline for processing eDNA bioinformatics workflows.

        :param input_path (str): The input path for raw sequences. This can be either a directory containing files or the path to a single file.
        :param output_path (str): The directory where output files will be saved.
        :param enabled_stages (list[str]): The process stages to be executed. The execution order will match the list order; each stage runs for all samples before the next one. Default:
         ["decompress", "merge", "cutprimer", "qfilter", "dereplicate", "denoise", "blast", "addlineage", "lca", "addhap"].
         "fqtofa" converts FASTQ to FASTA without filtering and can replace "qfilter".
        :param settings: Additional optional arguments to configure pipeline stages and runtime behavior. These include:
          Input File Suffix:
            - raw_suffix (str): File suffix for raw input sequences. Default: auto-detected from the starting stage.

          Merge Settings:
            - maxdiff (int): Maximum number of mismatches in the alignment. Default: 5.
            - pctid (int): Minimum %id of alignment. Default: 90.
            - min_merged_read_len (int): Discard merged reads shorter than this (-fastq_minmergelen). Default: None (no limit).
            - max_merged_read_len (int): Discard merged reads longer than this (-fastq_maxmergelen). Default: None (no limit).

          Cut Primer Settings:
            - rm_p_5 (str): Non-internal 5’ primer. Default: "GTCGGTAAAACTCGTGCCAGC" (MiFish-UF).
            - rm_p_3 (str): Non-internal 3’ primer. Default: "CAAACTGGGATTAGATACCCCACTATG" (reverse-complement MiFish-UR).
            - error_rate (float): The maximum rate of error could be tolerated. The actual error rate is computed as the number of errors in the match divided by the length of the matching part of the primer. Default: 0.15.
            - min_trimmed_read_len (int): Discard primer-trimmed reads shorter than this (written to cutprimer/discarded/<prefix>_too_short.fastq). Default: None (no limit; only empty reads are dropped).
            - max_trimmed_read_len (int): Discard primer-trimmed reads longer than this (written to cutprimer/discarded/<prefix>_too_long.fastq). Default: None (no limit).
            - revcomp (bool): If True, also search the reverse complement of each read for primer matches and output reads in the orientation where the primer was found. Default: True.

          Quality Filter Settings:
            - min_qscore (float): Keep reads whose expected number of errors per base (mean Phred error probability) is at most that of this quality score; converted to usearch -fastq_maxee_rate as 10^(-min_qscore/10), e.g. Q20 -> 0.01 (1.7 expected errors in a 170 bp read), Q30 -> 0.001. Default: 30.
            - max_qscore (int): Highest quality score in the reads (usearch -fastq_qmax); usearch stops on a higher score (93 is the highest FASTQ can store). Default: 70.
            - maxns (int): Discard reads with more Ns than this value. Default: None (no limit).

          Dereplicate Setting:
            - strand (bool): If True, treat a sequence and its reverse complement as duplicates during dereplication (-strand both). If False, only the forward strand is considered (-strand plus). Default: True.

          Denoise Settings:
            - minsize (int): Discard sequences with abundance that are smaller than this parameter. Default: 8.
            - alpha (int): Denoising sensitivity parameter. See UNOISE2 paper for definition. Default: 2.

          Cluster Setting:
            - cluster_minsize (int): Discard unique sequences with abundance smaller than this parameter before clustering. Default: 2.

          Cluster (ONT) Settings:
            - cluster_id (float): Identity for the initial vsearch clustering of raw reads. Default: 0.90.
            - min_cluster_size (int): Minimum reads in a cluster for it to be polished. Default: 10.
            - max_reads (int): Maximum reads per cluster used for polishing. Default: 200.
            - racon_rounds (int): Rounds of minimap2 + racon polishing. Default: 2.
            - medaka (bool): If True, polish with medaka after racon. Default: False.
            - medaka_model (str): Medaka model matching the basecaller (e.g. "r1041_e82_400bps_sup_v5.0.0"). Default: None (medaka's default).
            - consensus_id (float): Identity at which polished consensuses are merged into one OTU. Default: 0.99.
            - cluster_strand (bool): If True, cluster reads and merge consensuses on both strands. Default: False (reads are oriented by cutprimer).
            - keep_tmp (bool): If True, keep per-cluster intermediate files. Default: False.

          OTU Table Settings:
            - otutab_id (float): Minimum fractional identity for mapping sequences to OTUs/ZOTUs. Default: None (0.97 for denoise and cluster, 0.90 for cluster_ont, whose raw reads carry more errors).
            - otutab_strand (bool): If True, search both strands when mapping (-strand both). If False, only the plus strand is searched. Default: True.

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

          Add Haplotype Setting:
            - otu_dir (str): Path to the directory containing the OTU/ZOTU FASTA (<prefix>.fasta) and OTU table (<prefix>_otu_table.txt) used by addhap. Default: the output directory of the enabled denoise, cluster or cluster_ont stage (denoise if none is enabled).

          External Program Setting:
            - usearch_prog (str): Command to execute USEARCH for merge, qfilter, dereplicate, denoise, and cluster stages, and the OTU table of cluster_ont. Default: "usearch".
            - cutadapt_prog (str): Command to execute Cutadapt for cutprimer stage. Default: "cutadapt".
            - blast_prog (str): Command to execute BLAST for blast stage. Default: "blastn".
            - vsearch_prog, minimap2_prog, racon_prog, medaka_prog (str): Commands for the cluster_ont stage. Defaults: "vsearch", "minimap2", "racon", "medaka_consensus".

          Configuration Settings:
            - verbose (bool): Whether to enable verbose logging. Default: False.
            - dry (bool): If True, log every command ("DRY RUN: ...") without running it or writing output files; all stages are still set up for all samples. Default: False.
            - logger (Logger): Logger object for pipeline logging. Default: base_logger.logger.
            - n_cpu (int): Number of CPU cores to be used for processing. Default: 1.
            - memory (int): Maximum memory (in GB) allowed for processing. Default: 8.
        """
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"input path does not exist: {input_path}")
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
            "qfilter": qfilter.QFilterStage,
            "dereplicate": dereplicate.DereplicateStage,
            "denoise": denoise.DenoiseStage,
            "cluster": cluster.ClusterStage,
            "cluster_ont": cluster_ont.ClusterOntStage,
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
        self.qfilter_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["qfilter"].items()
        }
        self.dereplicate_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["dereplicate"].items()
        }
        self.denoise_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["denoise"].items()
        }
        self.cluster_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["cluster"].items()
        }
        self.cluster_ont_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["cluster_ont"].items()
        }
        self.otutab_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["otutab"].items()
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
        self.addhap_settings = {
            k: settings.get(k, v) for k, v in SETTINGS["addhap"].items()
        }
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
        if raw_suffix is not None:
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
        self.stage_in_dir = {}
        curr_dir = self.indir_path
        curr_suffix = self.stage_suffix["raw"]
        for stage in self.enabled_stages:
            if stage == "decompress":
                m = re.match(r".*(\.gz)$", curr_suffix)
                if m is not None:
                    self.stage_suffix["decompress"] = curr_suffix.replace(
                        ".gz", ""
                    )
                else:
                    self.config.logger.warning(
                        f"The in_suffix '{curr_suffix}' does not match the expected format (.gz) for the 'decompress' stage."
                        "Skipping the stage"
                    )
            elif stage in ["fqtofa", "qfilter"]:
                m = re.match(r".*\.(fastq|fq)$", curr_suffix)
                if m is not None:
                    self.stage_suffix[stage] = curr_suffix.replace(
                        m.group(1), "fasta"
                    )
                else:
                    self.config.logger.warning(
                        f"The in_suffix '{curr_suffix}' does not match the expected format (.fq/.fastq) for the '{stage}' stage."
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

            if stage in [
                "merge",
                "qfilter",
                "dereplicate",
                "denoise",
                "cluster",
                "cluster_ont",
            ]:
                stage_args["usearch_prog"] = self.prog_settings["usearch"]
            elif stage == "cutprimer":
                stage_args["cutadapt_prog"] = self.prog_settings["cutadapt"]
            elif stage in ["blast", "assigntaxa"]:
                stage_args["blast_prog"] = self.prog_settings["blast"]

            if stage == "addhap":
                stage_args["otu_dir"] = (
                    self.addhap_settings["otu_dir"]
                    or self.stage_dir[self.otu_stage()]
                )

            if stage in [
                "merge",
                "cutprimer",
                "qfilter",
                "dereplicate",
                "denoise",
                "blast",
                "addlineage",
                "lca",
                "cluster_ont",
            ]:
                stage_args.update(eval(f"self.{stage}_settings"))
            if stage == "cluster":
                stage_args["minsize"] = self.cluster_settings[
                    "cluster_minsize"
                ]
            if stage in ["denoise", "cluster", "cluster_ont"]:
                stage_args.update(
                    {
                        k: v
                        for k, v in self.otutab_settings.items()
                        if v is not None
                    }
                )
            if stage == "cluster_ont":
                for prog in [
                    "vsearch",
                    "minimap2",
                    "racon",
                    "medaka",
                ]:
                    stage_args[f"{prog}_prog"] = self.prog_settings[prog]

            self.stage_in_dir[stage] = curr_dir
            self.stages[stage] = self.stage_class[stage](**stage_args)

            curr_dir = self.stage_dir[stage]
            curr_suffix = self.stage_suffix[stage]

    def otu_stage(self):
        otu_stages = [
            s
            for s in self.enabled_stages
            if s in ["denoise", "cluster", "cluster_ont"]
        ]
        return otu_stages[-1] if otu_stages else "denoise"

    def run_stages(self, prefixes):
        for name, stage in self.stages.items():
            if not prefixes:
                self.config.logger.warning(
                    f"No samples left to run at stage: {name}"
                )
                return
            if getattr(stage, "pools_samples", False) is True:
                done = stage.run_all(prefixes)
                if self.config.verbose:
                    print()
                failed = [] if done else prefixes
            else:
                done = []
                for prefix in prefixes:
                    if self.run_stage(name, stage, prefix):
                        done.append(prefix)
                failed = [p for p in prefixes if p not in done]
            for prefix in failed:
                self.config.logger.error(
                    f"Error: process errors at stage: {name} ({prefix})\n"
                )
            prefixes = done

    def run_stage(self, name, stage, prefix):
        base_logger.logger.info(f"Sample ID: {prefix}")
        try:
            stage.setup(prefix)
        except (FileNotFoundError, ValueError) as e:
            self.config.logger.error(f"FAIL: setup of stage {name}: {e}")
            stage.runners = []
            return False
        is_complete = stage.run()
        if self.config.verbose:
            print()
        return is_complete

    def run_stages_files(self):
        suffix = self.stage_suffix["raw"]
        side_suffixes = [
            s for s in SIDE_FILE_SUFFIXES if s.endswith(suffix) and s != suffix
        ]
        prefixes = sorted(
            file[: -len(suffix)]
            for file in os.listdir(self.indir_path)
            if file.endswith(suffix)
            and not any(file.endswith(s) for s in side_suffixes)
        )

        if len(prefixes) == 0:
            self.config.logger.warning(
                f"No input files found with suffix '{suffix}' in directory '{self.indir_path}'"
            )
            return

        self.run_stages(prefixes)

    def run_stages_one_file(self, input_path):
        prefix = os.path.basename(input_path).replace(
            self.stage_suffix["raw"], ""
        )
        self.run_stages([prefix])

    def close_file_handler(self):
        base_logger.close_file_handler(self.config.logger)
