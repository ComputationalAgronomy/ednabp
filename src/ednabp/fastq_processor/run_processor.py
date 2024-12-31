import os

from ..analysis_toolkit.runner_build import base_logger
from .step_build import stage_config
from .step_exec import (decompress, merge, cut_primer, fq_to_fa, dereplicate, denoise, assign_taxa,)

class FastqProcessor:

    def __init__(self,
                 input_path: str,
                 output_path: str,
                 **settings
                 ):
        '''
        A pipeline for processing eDNA bioinformatics workflows.
        
        :param input_path (str): The input path for raw sequences. This can be either a directory containing files or the path to a single file.
        :param output_path (str): The directory where output files will be saved.
        :param enabled_stages (list[str]): The process stages to be executed. The execution order will match the list order. Default:
         ["decompress", "merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "assigntaxa"] (all stages will be run).
        :param settings: Additional optional arguments to configure pipeline stages and runtime behavior. These include:
          Directory Names:
            - decompress_dir_name (str): The name of the subdirectory for the decompression stage. Default: "decompress".
            - merge_dir_name (str): The name of the subdirectory for the merge stage. Default: "merge".
            - cutprimer_dir_name (str): The name of the subdirectory for the cut-primer stage. Default: "cutprimer".
            - fqtofa_dir_name (str): The name of the subdirectory for converting FastQ to FastA. Default: "fqtofa".
            - dereplicate_dir_nam (str): The name of the subdirectory for the dereplication stage. Default: "dereplicate".
            - denoise_dir_name (str): The name of the subdirectory for the denoising stage. Default: "denoise".
            - assigntaxa_dir_name (str): The name of the subdirectory for taxonomic assignment. Default: "assigntaxa".

          File Suffixes:
            - raw_suffix (str): File suffix for raw input sequences. Default: "_R1.fastq.gz".
            - decompress_suffix (str): File suffix for sequences after decompression. Default: "_R1.fastq".
            - merge_suffix (str): File suffix for merged sequences. Default: "_merged.fastq".
            - cutprimer_suffix (str): File suffix for trimmed sequences after primer removal. Default: "_trimmed.fastq".
            - dereplicate_suffix (str): File suffix for unique sequences after dereplication. Default: "_uniqs.fasta".
            - denoise_suffix (str): File suffix for denoised sequences (ZOTUs). Default: "_zotus.fasta".
            - assigntaxa_suffix (str): File suffix for taxonomic assignment results. Default: "_taxa.csv".

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

          Assign Taxa Settings:
            - db_path (str): Path to the taxonomic database. Default: None.
            - lineage_path (str): Path to the taxonomic lineage file. Default: None.
            - evalue (float): Expectation value (E) threshold for saving hits. Default: 0.00001.
            - qcov_hsp_perc (int): The %threshold of the query sequence that has to form an alignment against the reference to be retained. Default: 90.
            - perc_identity (int): Minimum percentage identity required for taxonomic assignment. Default: 90.
            - specifiers (str): Output format specifiers for BLAST results. Default: 
              "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".

          Configuration Settings:
            - verbose (bool): Whether to enable verbose logging. Default: True.
            - dry (bool): If True, perform a dry run without executing commands. Default: False.
            - logger (Logger): Logger object for pipeline logging. Default: base_logger.logger.
            - n_cpu (int): Number of CPU cores to be used for processing. Default: 1.
            - memory (int): Maximum memory (in GB) allowed for processing. Default: 8.
        '''
        assert os.path.exists(input_path), f"Error: input path does not exist: {input_path}"
        os.makedirs(output_path, exist_ok=True)
        input_is_dir = True if os.path.isdir(input_path) else False

        self.indir_path = input_path if input_is_dir else os.path.dirname(input_path)
        self.outdir_path = output_path

        self.add_default_settings(settings)

        self.add_config()

        self.setup_stages()

        if input_is_dir:
            self.run_stages_files()
        else:
            self.run_stages_one_file(input_path)

    def add_default_settings(self, settings):
        self.stage_class = {
            "decompress": decompress.DecompressStage,
            "merge": merge.MergeStage,
            "cutprimer": cut_primer.CutPrimerStage,
            "fqtofa": fq_to_fa.FqToFaStage,
            "dereplicate": dereplicate.DereplicateStage,
            "denoise": denoise.DenoiseStage,
            "assigntaxa": assign_taxa.AssignTaxaStage,
        }

        DEFAULT_SETTINGS = {
            'enabled_stages': ["decompress", "merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "assigntaxa"],
            'stage_dir_name':{
                "decompress": "decompress",
                "merge": "merge",
                "cutprimer": "cutprimer",
                "fqtofa": "fqtofa",
                "dereplicate": "dereplicate",
                "denoise": "denoise",
                "assigntaxa": "assigntaxa",
            },
            'suffix': {
                "raw": "_R1.fastq.gz",
                "decompress": "_R1.fastq",
                "merge": "_merged.fastq",
                "cutprimer": "_trimmed.fastq",
                "dereplicate": "_uniqs.fasta",
                "denoise": "_zotus.fasta",
                "assigntaxa": "_taxa.csv",
            },
            'merge': {
                "maxdiff": 5,
                "pctid": 90,
            },
            'cutprimer': {
                "rm_p_5": "GTCGGTAAAACTCGTGCCAGC",
                "rm_p_3": "CAAACTGGGATTAGATACCCCACTATG",
                "error_rate": 0.15,
                "min_read_len": 204,
                "max_read_len": 254,
            },
            'denoise': {
                "minsize": 8,
                "alpha": 2,
            },
            'assigntaxa': {
                "db_path": None,
                "lineage_path": None,
                "evalue": 0.00001,
                "qcov_hsp_perc": 90,
                "perc_identity": 90,
                "specifiers": "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            },
            'config': {
                "verbose": True,
                "dry": False,
                "logger": base_logger.logger,
                "n_cpu": 1,
                "memory": 8
            }
        }

        self.enabled_stages = settings.get("enabled_stages", DEFAULT_SETTINGS['enabled_stages'])
        self.stage_dir_name = {k: settings.get(f"{k}_dir_name", v) for k, v in DEFAULT_SETTINGS['stage_dir_name'].items()}
        self.stage_dir = {k: os.path.join(self.outdir_path, v) for k, v in self.stage_dir_name.items()}
        self.stage_suffix = {k: settings.get(f"{k}_suffix", v) for k, v in DEFAULT_SETTINGS['suffix'].items()}
        self.merge_settings = {k: settings.get(k, v) for k, v in DEFAULT_SETTINGS['merge'].items()}
        self.cutprimer_settings = {k: settings.get(k, v) for k, v in DEFAULT_SETTINGS['cutprimer'].items()}
        self.denoise_settings = {k: settings.get(k, v) for k, v in DEFAULT_SETTINGS['denoise'].items()}
        self.assigntaxa_settings = {k: settings.get(k, v) for k, v in DEFAULT_SETTINGS['assigntaxa'].items()}
        self.config_settings = {k: settings.get(k, v) for k, v in DEFAULT_SETTINGS['config'].items()}

    def add_config(self):
        fp_fh = base_logger._get_file_handler(os.path.join(self.outdir_path, "stages.log"))
        self.config_settings["logger"].addHandler(fp_fh)
        self.config = stage_config.StageConfig(settings = self.config_settings)

    def setup_stages(self):
        self.stages = dict()
        curr_dir = self.indir_path
        curr_suffix = self.stage_suffix["raw"]
        for stage in self.enabled_stages:
            if stage == "fqtofa":
                self.stage_suffix["fqtofa"] = curr_suffix.replace("fastq", "fasta")
            stage_args = {
                "config": self.config,
                "in_dir": curr_dir,
                "out_dir": self.stage_dir[stage],
                "in_suffix": curr_suffix,
                "out_suffix": self.stage_suffix[stage],
            }
            if stage in ["merge", "cutprimer", "denoise", "assigntaxa"]:
                stage_args.update(eval(f"self.{stage}_settings"))

            self.stages[stage] = self.stage_class[stage](**stage_args)

            curr_dir = self.stage_dir[stage]
            curr_suffix = self.stage_suffix[stage]

    def run_one_file(self, prefix):
        print(f"Sample ID: {prefix}")
        for k, s in self.stages.items():
            s.setup(prefix)
            is_complete = s.run()
            if not is_complete:
                print(f"Error: process errors at stage: {k}\n")
                break
            print()

    def run_stages_files(self):
        suffix = self.stage_suffix["raw"]
        files = os.listdir(self.indir_path)
        prefixes = [file.replace(suffix, "") for file in files if file.endswith(suffix)]
        for prefix in prefixes:
            self.run_one_file(prefix)

    def run_stages_one_file(self, input_path):
        prefix = os.path.basename(input_path).replace(self.stage_suffix["raw"], "")
        self.run_one_file(prefix, self.stages)