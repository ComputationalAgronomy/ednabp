import os

from ..analysis_toolkit.runner_build import base_logger
from .step_build import stage_config
from .step_exec import (decompress, merge, cut_primer, fq_to_fa, dereplicate, denoise, assign_taxa,)

class FastqProcessor:

    @staticmethod
    def get_prefix_with_suffix(in_dir, suffix):
        files = os.listdir(in_dir)
        prefix = [file.replace(suffix, "") for file in files if file.endswith(suffix)]
        return prefix

    @staticmethod
    def run_single_data(prefix, stages):
        print(f"Sample ID: {prefix}")
        for k, s in stages.items():
            s.setup(prefix)
            is_complete = s.run()
            if not is_complete:
                print(f"Error: process errors at stage: {k}\n")
                break
            print()

    def __init__(self,
                 input_path: str,
                 output_path: str,
                 enabled_stages=["decompress", "merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "assigntaxa"],
                 **settings
                 ):
        assert os.path.exists(input_path), f"Error: input path does not exist: {input_path}"
        os.makedirs(output_path, exist_ok=True)
        input_is_dir = True if os.path.isdir(input_path) else False

        self.indir_path = input_path if input_is_dir else os.path.dirname(input_path)
        self.outdir_path = output_path

        self.add_default_settings(settings)

        fp_fh = base_logger._get_file_handler(os.path.join(output_path, "stages.log"))
        self.config_settings["logger"].addHandler(fp_fh)
        self.config = stage_config.StageConfig(settings = self.config_settings)

        self.setup_stages(enabled_stages)

        if input_is_dir:
            self.data_prefix = FastqProcessor.get_prefix_with_suffix(self.indir_path, self.stage_suffix["raw"])
            for prefix in self.data_prefix:
                FastqProcessor.run_single_data(prefix, self.stages)
        else:
            prefix = os.path.basename(input_path).replace(self.stage_suffix["raw"], "")
            self.run_single_data(prefix, self.stages)

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

        self.stage_dir_name = {k: settings.pop(f"{k}_dir_name", v) for k, v in DEFAULT_SETTINGS['stage_dir_name'].items()}
        self.stage_dir = {k: os.path.join(self.outdir_path, v) for k, v in self.stage_dir_name.items()}
        self.stage_suffix = {k: settings.pop(f"{k}_suffix", v) for k, v in DEFAULT_SETTINGS['suffix'].items()}
        self.merge_settings = {k: settings.pop(k, v) for k, v in DEFAULT_SETTINGS['merge'].items()}
        self.cutprimer_settings = {k: settings.pop(k, v) for k, v in DEFAULT_SETTINGS['cutprimer'].items()}
        self.denoise_settings = {k: settings.pop(k, v) for k, v in DEFAULT_SETTINGS['denoise'].items()}
        self.assigntaxa_settings = {k: settings.pop(k, v) for k, v in DEFAULT_SETTINGS['assigntaxa'].items()}
        self.config_settings = {k: settings.pop(k, v) for k, v in DEFAULT_SETTINGS['config'].items()}

    def setup_stages(self,
                     enabled_stages,
                     ):
        self.stages = dict()
        curr_dir = self.indir_path
        curr_suffix = self.stage_suffix["raw"]
        for stage in enabled_stages:
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

def main():
    FastqProcessor(input_path=".\\stage_test\\fastq",
                   output_path=".\\stage_test",
                   db_path=".\\data\\database\\MiFish",
                   lineage_path=".\\data\\database\\lineage.csv",
                   n_cpu=20,
                   )

if __name__ == "__main__":
    main()