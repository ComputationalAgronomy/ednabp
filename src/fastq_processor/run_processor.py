import os

from analysis_toolkit.runner_build import base_logger
from fastq_processor.step_build import stage_config
from fastq_processor.step_exec import (decompress, merge, cut_primer, fq_to_fa, 
                                       dereplicate, denoise, assign_taxa,)

class FastqProcessor:

    @staticmethod
    def get_prefix_with_suffix(in_dir, suffix):
        files = os.listdir(in_dir)
        prefix = [file.replace(suffix, "") for file in files if file.endswith(suffix)]
        return prefix

    @staticmethod
    def run_each_data(prefix, stages):
        print(f"Sample ID: {prefix}")
        for k, s in stages.items():
            s.setup(prefix)
            is_complete = s.run()
            if not is_complete:
                print(f"Error: process errors at stage: {k}\n")
                break
            print()

    def __init__(self,
        stages_parent_dir: str,
        raw_dir_name: str,
        db_path: str,
        lineage_path: str,
        enabled_stages=["decompress", "merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "assigntaxa"],
        decompress_dir_name: str = "decompress",
        merge_dir_name: str = "merge",
        cutprimer_dir_name: str = "cut_primer",
        fqtofa_dir_name: str = "fq_to_fa",
        derep_dir_name: str = "dereplicate",
        denoise_dir_name: str = "denoise",
        blast_dir_name: str = "blast",
        raw_suffix: str = "R1.fastq.gz",
        decompress_suffix: str = "R1.fastq",
        merge_suffix: str = "merge.fastq",
        cutprimer_suffix: str = "cut.fastq",
        derep_suffix: str = "uniq.fasta",
        denoise_suffix: str = "denoise.fasta",
        blast_suffix: str = "blast.csv",
        maxdiff: int = 5,
        pctid: int = 90,
        rm_p_5: str = "GTCGGTAAAACTCGTGCCAGC",
        rm_p_3: str = "CAAACTGGGATTAGATACCCCACTATG",
        error_rate: float = 0.15,
        min_read_len: int = 204,
        max_read_len: int = 254,
        minsize: int = 8,
        alpha: int = 2,
        evalue: float = 0.00001,
        qcov_hsp_perc: int = 90,
        perc_identity: int = 90,
        specifiers: str = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        verbose: bool=True,
        n_cpu: int = 1,
        memory: int = 8,
        ):
        fp_fh = base_logger._get_file_handler(os.path.join(stages_parent_dir, "stages.log"))
        base_logger.logger.addHandler(fp_fh)
        self.config = stage_config.StageConfig(verbose=verbose, logger=base_logger.logger, n_cpu=n_cpu, memory=memory) # TODO(SW): Expand this class, so you only need to pass one parameters.
        self.parent_dir = stages_parent_dir
        self.input_dir = os.path.join(stages_parent_dir, raw_dir_name)
        self.setup_stages(
            enabled_stages=enabled_stages,
            stages_parent_dir=stages_parent_dir,
            raw_dir_name=raw_dir_name,
            decompress_dir_name=decompress_dir_name,
            merge_dir_name=merge_dir_name,
            cutprimer_dir_name=cutprimer_dir_name,
            fqtofa_dir_name=fqtofa_dir_name,
            derep_dir_name=derep_dir_name,
            denoise_dir_name=denoise_dir_name,
            blast_dir_name=blast_dir_name,
            raw_suffix=raw_suffix,
            decompress_suffix=decompress_suffix,
            merge_suffix=merge_suffix,
            cutprimer_suffix=cutprimer_suffix,
            derep_suffix=derep_suffix,
            denoise_suffix=denoise_suffix,
            blast_suffix=blast_suffix,
            db_path=db_path,
            lineage_path=lineage_path,
            maxdiff=maxdiff,
            pctid=pctid,
            rm_p_5=rm_p_5,
            rm_p_3=rm_p_3,
            error_rate=error_rate,
            min_read_len=min_read_len,
            max_read_len=max_read_len,
            minsize=minsize,
            alpha=alpha,
            evalue=evalue,
            qcov_hsp_perc=qcov_hsp_perc,
            perc_identity=perc_identity,
            specifiers=specifiers
        )
        self.data_prefix = FastqProcessor.get_prefix_with_suffix(self.input_dir, f"_{raw_suffix}")

        for prefix in self.data_prefix:
            FastqProcessor.run_each_data(prefix, self.stages)

    def setup_stages(self,
            enabled_stages,
            stages_parent_dir,
            raw_dir_name,
            decompress_dir_name,
            merge_dir_name,
            cutprimer_dir_name,
            fqtofa_dir_name,
            derep_dir_name,
            denoise_dir_name,
            blast_dir_name,
            raw_suffix,
            decompress_suffix,
            merge_suffix,
            cutprimer_suffix,
            derep_suffix,
            denoise_suffix,
            blast_suffix,
            db_path,
            lineage_path,
            maxdiff,
            pctid,
            rm_p_5,
            rm_p_3,
            error_rate,
            min_read_len,
            max_read_len,
            minsize,
            alpha,
            evalue,
            qcov_hsp_perc,
            perc_identity,
            specifiers
        ):
        raw_dir = os.path.join(stages_parent_dir, raw_dir_name)
        decompress_dir = os.path.join(stages_parent_dir, decompress_dir_name)
        merge_dir = os.path.join(stages_parent_dir, merge_dir_name)
        cutprimer_dir = os.path.join(stages_parent_dir, cutprimer_dir_name)
        fqtofa_dir = os.path.join(stages_parent_dir, fqtofa_dir_name)
        derep_dir = os.path.join(stages_parent_dir, derep_dir_name)
        denoise_dir = os.path.join(stages_parent_dir, denoise_dir_name)
        blast_dir = os.path.join(stages_parent_dir, blast_dir_name)

        self.stages = dict()
        curr_dir = raw_dir
        curr_suffix = raw_suffix
        if "decompress" in enabled_stages:
            self.stages["decompress"] = decompress.DecompressStage(self.config, in_dir=curr_dir, out_dir=decompress_dir,
                                                                   in_suffix=curr_suffix, out_suffix=decompress_suffix)
            curr_dir = decompress_dir
            curr_suffix = decompress_suffix
        if "merge" in enabled_stages:
            self.stages["merge"] = merge.MergeStage(self.config, in_dir=curr_dir, out_dir=merge_dir,
                                                    in_suffix=curr_suffix, out_suffix=merge_suffix,
                                                    maxdiff=maxdiff, pctid=pctid)
            curr_dir = merge_dir
            curr_suffix = merge_suffix
        if "cutprimer" in enabled_stages:
            self.stages["cutprimer"] = cut_primer.CutPrimerStage(self.config, in_dir=curr_dir, out_dir=cutprimer_dir,
                                                                 in_suffix=curr_suffix, out_suffix=cutprimer_suffix,
                                                                 rm_p_5=rm_p_5, rm_p_3=rm_p_3, error_rate=error_rate,
                                                                 min_read_len=min_read_len, max_read_len=max_read_len)
            curr_dir = cutprimer_dir
            curr_suffix = cutprimer_suffix
        if "fqtofa" in enabled_stages:
            out_suffix = curr_suffix.replace("fastq", "fasta")
            self.stages["fqtofa"] = fq_to_fa.FqToFaStage(self.config, in_dir=curr_dir, out_dir=fqtofa_dir,
                                                         in_suffix=curr_suffix, out_suffix=out_suffix)
            curr_dir = fqtofa_dir
            curr_suffix = out_suffix
        if "dereplicate" in enabled_stages:
            self.stages["dereplicate"] = dereplicate.DereplicateStage(self.config, in_dir=curr_dir, out_dir=derep_dir,
                                                                      in_suffix=curr_suffix, out_suffix=derep_suffix)
            curr_dir = derep_dir
            curr_suffix = derep_suffix
        if "denoise" in enabled_stages:
            self.stages["denoise"] = denoise.DenoiseStage(self.config, in_dir=curr_dir, out_dir=denoise_dir,
                                                          in_suffix=curr_suffix, out_suffix=denoise_suffix,
                                                          minsize=minsize, alpha=alpha)
            curr_dir = denoise_dir
            curr_suffix = denoise_suffix
        if "assigntaxa" in enabled_stages:
            self.stages["assigntaxa"] = assign_taxa.AssignTaxaStage(self.config, in_dir=curr_dir, out_dir=blast_dir,
                                                                    in_suffix=curr_suffix, out_suffix=blast_suffix,
                                                                    db_path=db_path, lineage_path=lineage_path,
                                                                    evalue=evalue, qcov_hsp_perc=qcov_hsp_perc,
                                                                    perc_identity=perc_identity, specifiers=specifiers)

def main():
    FastqProcessor(
        stages_parent_dir="stage_test",
        fastq_dir_name="fastq",
        db_path="data\\database\\MiFish",
        lineage_path="data\\database\\lineage.csv",
        n_cpu=20,
    )

if __name__ == "__main__":
    main()