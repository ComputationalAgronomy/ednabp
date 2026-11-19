#!/usr/bin/env python

import argparse

from ..common.default_settings import SETTINGS


def main(argv=None):
    argument_groups = {
        "Required Input & Output Path": [
            {
                "args": ["-i", "--input-path"],
                "kwargs": {
                    "required": True,
                    "help": "Input FASTA sequence file path or folder path.",
                },
            },
            {
                "args": ["-o", "--output-path"],
                "kwargs": {
                    "required": True,
                    "help": "Directory where output files will be saved.",
                },
            },
        ],
        "Stages Executed": [
            {
                "args": ["-e", "--enabled-stages"],
                "kwargs": {
                    "nargs": "+",
                    "type": str,
                    "default": SETTINGS["enabled_stages"],
                    "help": "Process stages to be executed in order; each stage runs for all samples before the next. Add 'pool' after dereplicate (pool before OTU generation) or after denoise/cluster/cluster_ont (pool the samples' OTUs); the stages after it run once on the pooled sample (default [run all stages]: %(default)s).",
                },
            },
        ],
        "Input File Suffix": [
            {
                "args": ["--raw-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["raw"],
                    "help": "File suffix for raw input files (default: auto-detected from starting stage).",
                },
            },
        ],
        "Paired-end Merge Settings": [
            {
                "args": ["--maxdiff"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["merge"]["maxdiff"],
                    "help": "Maximum number of mismatches in alignment during merging (default: %(default)s).",
                },
            },
            {
                "args": ["--pctid"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["merge"]["pctid"],
                    "help": "Minimum percentage identity required for merging reads (default: %(default)s).",
                },
            },
            {
                "args": ["--min-merged-read-len"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["merge"]["min_merged_read_len"],
                    "help": "Discard merged reads shorter than this length (default: no limit).",
                },
            },
            {
                "args": ["--max-merged-read-len"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["merge"]["max_merged_read_len"],
                    "help": "Discard merged reads longer than this length (default: no limit).",
                },
            },
        ],
        "Primer Trimming Settings": [
            {
                "args": ["--rm-p-5"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["cutprimer"]["rm_p_5"],
                    "help": "Non-internal 5’ primer (default [MiFish-UF]: %(default)s).",
                },
            },
            {
                "args": ["--rm-p-3"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["cutprimer"]["rm_p_3"],
                    "help": "Non-internal 3’ primer (default [rev-com MiFish-UR]: %(default)s).",
                },
            },
            {
                "args": ["--error-rate"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["cutprimer"]["error_rate"],
                    "help": "Maximum allowable error rate for primer matching (default: %(default)s).",
                },
            },
            {
                "args": ["--min-trimmed-read-len"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cutprimer"]["min_trimmed_read_len"],
                    "help": "Discard primer-trimmed reads shorter than this length (default: no limit).",
                },
            },
            {
                "args": ["--max-trimmed-read-len"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cutprimer"]["max_trimmed_read_len"],
                    "help": "Discard primer-trimmed reads longer than this length (default: no limit).",
                },
            },
            {
                "args": ["--revcomp"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["cutprimer"]["revcomp"],
                    "help": "Also search the reverse complement of each read for primer matches and output reads in the matched orientation (default: %(default)s).",
                },
            },
        ],
        "Quality Filter Settings": [
            {
                "args": ["--min-qscore"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["qfilter"]["min_qscore"],
                    "help": "Keep reads whose expected errors per base are at most those of this quality score, e.g. 20 = 1 error per 100 bases (default: %(default)s).",
                },
            },
            {
                "args": ["--max-qscore"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["qfilter"]["max_qscore"],
                    "help": "Highest quality score in the reads; any value at least that high works (default: %(default)s).",
                },
            },
            {
                "args": ["--maxns"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["qfilter"]["maxns"],
                    "help": "Discard reads with more Ns than this value (default: no limit).",
                },
            },
        ],
        "Dereplicate Setting": [
            {
                "args": ["--strand"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["dereplicate"]["strand"],
                    "help": "Treat a sequence and its reverse complement as duplicates during dereplication (default: %(default)s).",
                },
            },
        ],
        "Denoise Settings": [
            {
                "args": ["--minsize"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["denoise"]["minsize"],
                    "help": "Minimum abundance of sequences to retain during denoising (default: %(default)s).",
                },
            },
            {
                "args": ["--alpha"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["denoise"]["alpha"],
                    "help": "Denoising sensitivity parameter (default: %(default)s).",
                },
            },
        ],
        "Cluster Setting": [
            {
                "args": ["--cluster-minsize"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cluster"]["cluster_minsize"],
                    "help": "Minimum abundance of unique sequences to retain before UPARSE OTU clustering; used by the cluster stage, an alternative to denoise (default: %(default)s).",
                },
            },
        ],
        "Cluster (ONT) Settings": [
            {
                "args": ["--cluster-id"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["cluster_ont"]["cluster_id"],
                    "help": "Identity for the initial vsearch clustering of raw ONT reads in the cluster_ont stage (default: %(default)s).",
                },
            },
            {
                "args": ["--min-cluster-size"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cluster_ont"]["min_cluster_size"],
                    "help": "Minimum reads in a cluster for it to be polished (default: %(default)s).",
                },
            },
            {
                "args": ["--max-reads"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cluster_ont"]["max_reads"],
                    "help": "Maximum reads per cluster used for polishing (default: %(default)s).",
                },
            },
            {
                "args": ["--racon-rounds"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cluster_ont"]["racon_rounds"],
                    "help": "Rounds of minimap2 + racon polishing (default: %(default)s).",
                },
            },
            {
                "args": ["--medaka"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["cluster_ont"]["medaka"],
                    "help": "Polish with medaka after racon (default: %(default)s).",
                },
            },
            {
                "args": ["--medaka-model"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["cluster_ont"]["medaka_model"],
                    "help": "Medaka model matching the basecaller, e.g. r1041_e82_400bps_sup_v5.0.0 (default: %(default)s).",
                },
            },
            {
                "args": ["--consensus-id"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["cluster_ont"]["consensus_id"],
                    "help": "Identity at which polished consensus sequences are merged into one OTU (default: %(default)s).",
                },
            },
            {
                "args": ["--cluster-strand"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["cluster_ont"]["cluster_strand"],
                    "help": "Cluster reads and merge consensus sequences on both strands in the cluster_ont stage; the abundance table uses --otutab-strand (default: %(default)s).",
                },
            },
            {
                "args": ["--keep-tmp"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["cluster_ont"]["keep_tmp"],
                    "help": "Keep per-cluster intermediate files (default: %(default)s).",
                },
            },
        ],
        "Pool Setting": [
            {
                "args": ["--pool-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["pool"]["pool_name"],
                    "help": "Name of the pooled sample written by the pool stage (default: %(default)s).",
                },
            },
        ],
        "OTU Table Settings": [
            {
                "args": ["--otutab-id"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["otutab"]["otutab_id"],
                    "help": "Minimum fractional identity for mapping sequences to OTUs/ZOTUs when building the abundance table (default: 0.97 for denoise and cluster, 0.90 for cluster_ont).",
                },
            },
            {
                "args": ["--otutab-strand"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["otutab"]["otutab_strand"],
                    "help": "Search both strands when mapping sequences to OTUs/ZOTUs for the abundance table in the denoise, cluster and cluster_ont stages (default: %(default)s).",
                },
            },
        ],
        "BLAST Settings": [
            {
                "args": ["--evalue"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["blast"]["evalue"],
                    "help": "Expectation value (E) threshold for saving hits (default: %(default)s).",
                },
            },
            {
                "args": ["--qcov-hsp-perc"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["blast"]["qcov_hsp_perc"],
                    "help": "Percentage of query coverage for high-scoring segment pairs (default: %(default)s).",
                },
            },
            {
                "args": ["--perc-identity"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["blast"]["perc_identity"],
                    "help": "Minimum Percentage identity for taxonomic assignment (default: %(default)s).",
                },
            },
            {
                "args": ["--maxhitnum"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["blast"]["maxhitnum"],
                    "help": "Maximum number of hits to keep per query sequence (default: %(default)s).",
                },
            },
            {
                "args": ["--specifiers"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["blast"]["specifiers"],
                    "help": "Custom format specifiers for BLAST results. (default: %(default)s).",
                },
            },
            {
                "args": ["-bdb", "--blast-db"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["blast"]["blast_db"],
                    "help": "Path to the BLAST database (default [NCBI remote nucleotide sequence database]: %(default)s).",
                },
            },
            {
                "args": ["--task"],
                "kwargs": {
                    "type": str,
                    "choices": [
                        "megablast",
                        "dc-megablast",
                        "blastn",
                        "blastn-short",
                    ],
                    "default": SETTINGS["blast"]["task"],
                    "help": "BLAST search task: megablast for highly similar sequences, dc-megablast for more dissimilar ones, blastn for sensitive search, blastn-short for queries shorter than 50 bp (default: %(default)s).",
                },
            },
            {
                "args": ["--taxids"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["blast"]["taxids"],
                    "help": "Comma-separated NCBI taxids to limit a local BLAST search to, including their child taxa, e.g. 2759 (eukaryotes) or 33208,33090 (animals, plants); needs taxonomy4blast.sqlite3 next to the database, ignored for remote nt (default: all taxa).",
                },
            },
        ],
        "Add Lineage Settings": [
            {
                "args": ["-lndb", "--lineage-db"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["addlineage"]["lineage_db"],
                    "help": "Path to the lineage database (default [NCBI remote nucleotide lineage database]: %(default)s).",
                },
            },
            {
                "args": ["-email", "--entrez-email"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["addlineage"]["entrez_email"],
                    "help": "The email used by NCBI to contact you in case of excessive usage or issues. (default : %(default)s).",
                },
            },
        ],
        "LCA Settings": [
            {
                "args": ["--tol-pct"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["lca"]["tol_pct"],
                    "help": "Percentage of the top bitscore used as the inclusion threshold before LCA. Hits with bitscore >= top * (1 - tol_pct / 100) are included in the consensus (default: %(default)s).",
                },
            },
            {
                "args": ["--score-column"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["lca"]["score_column"],
                    "help": "Column used for score-based hit filtering (default: %(default)s).",
                },
            },
            {
                "args": ["--qseqid-column"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["lca"]["qseqid_column"],
                    "help": "Column used to group hits by query sequence (default: %(default)s).",
                },
            },
        ],
        "Add Haplotype Setting": [
            {
                "args": ["--otu-dir"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["addhap"]["otu_dir"],
                    "help": "Directory with the OTU/ZOTU FASTA and OTU table used by addhap (default: output directory of the enabled denoise, cluster or cluster_ont stage).",
                },
            },
        ],
        "External Program Settings": [
            {
                "args": ["--usearch-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["usearch"],
                    "help": "Command to execute USEARCH for merge, qfilter, dereplicate, denoise, and cluster stages, and the OTU table of cluster_ont (default: %(default)s).",
                },
            },
            {
                "args": ["--cutadapt-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["cutadapt"],
                    "help": "Command to execute Cutadapt for cutprimer stage (default: %(default)s).",
                },
            },
            {
                "args": ["--blast-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["blast"],
                    "help": "Command to execute BLAST for blast stage (default: %(default)s).",
                },
            },
            {
                "args": ["--vsearch-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["vsearch"],
                    "help": "Command to execute VSEARCH for cluster_ont stage (default: %(default)s).",
                },
            },
            {
                "args": ["--minimap2-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["minimap2"],
                    "help": "Command to execute minimap2 for cluster_ont stage (default: %(default)s).",
                },
            },
            {
                "args": ["--racon-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["racon"],
                    "help": "Command to execute Racon for cluster_ont stage (default: %(default)s).",
                },
            },
            {
                "args": ["--medaka-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["medaka"],
                    "help": "Command to execute medaka_consensus for cluster_ont stage (default: %(default)s).",
                },
            },
        ],
        "Configuration Settings": [
            {
                "args": ["--verbose"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["config_basic"]["verbose"],
                    "help": "Enable detailed logging output (default: %(default)s).",
                },
            },
            {
                "args": ["--dry"],
                "kwargs": {
                    "action": argparse.BooleanOptionalAction,
                    "default": SETTINGS["config_basic"]["dry"],
                    "help": "Log every command ('DRY RUN: ...') without running it or writing output files, to check a run before starting it (default: %(default)s).",
                },
            },
            {
                "args": ["--n-cpu"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["config_machine"]["n_cpu"],
                    "help": "Number of CPU cores to use for processing (default: %(default)s).",
                },
            },
        ],
    }

    parser = argparse.ArgumentParser(
        formatter_class=lambda prog: MyFormatter(prog, max_help_position=32),
        argument_default=argparse.SUPPRESS,
        description="A pipeline for processing eDNA bioinformatics workflows.",
    )

    for group_name, args_list in argument_groups.items():
        group = parser.add_argument_group(group_name)
        for arg_data in args_list:
            group.add_argument(*arg_data["args"], **arg_data["kwargs"])

    options = vars(parser.parse_args(argv))

    from ..bp import BioPipeline

    try:
        BioPipeline(**options)
    except FileNotFoundError as e:
        parser.error(str(e))


class MyFormatter(argparse.HelpFormatter):
    """
    for matt wilkie on SO
    https://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            (metavar,) = self._metavar_formatter(action, default)(1)
            return metavar
        else:
            parts = []
            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)
            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                args_string = self._format_args(action, "")
                for option_string in action.option_strings:
                    parts.append(option_string)
                return "{} {}".format(", ".join(parts), args_string)
            return ", ".join(parts)


if __name__ == "__main__":
    main()
