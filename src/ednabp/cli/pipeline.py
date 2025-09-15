#!/usr/bin/env python

import argparse

from ..common.default_settings import SETTINGS


def main():
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
                    "help": "Process stages to be executed in order (default [run all stages]: %(default)s).",
                },
            }
        ],
        "Directory Names": [
            {
                "args": ["--decompress-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["decompress"],
                    "help": "Name of the subdirectory for the decompression stage (default: %(default)s).",
                },
            },
            {
                "args": ["--merge-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["merge"],
                    "help": "Name of the subdirectory for the merge stage (default: %(default)s).",
                },
            },
            {
                "args": ["--cutprimer-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["cutprimer"],
                    "help": "Name of the subdirectory for the cut-primer stage (default: %(default)s).",
                },
            },
            {
                "args": ["--fqtofa-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["fqtofa"],
                    "help": "Name of the subdirectory for FastQ to FastA conversion (default: %(default)s).",
                },
            },
            {
                "args": ["--dereplicate-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["dereplicate"],
                    "help": "Name of the subdirectory for the dereplication stage (default: %(default)s).",
                },
            },
            {
                "args": ["--denoise-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["denoise"],
                    "help": "Name of the subdirectory for the denoising stage (default: %(default)s).",
                },
            },
            {
                "args": ["--assigntaxa-dir-name"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["dir_name"]["assigntaxa"],
                    "help": "Name of the subdirectory for taxonomic assignment (default: %(default)s).",
                },
            },
        ],
        "File Suffixes": [
            {
                "args": ["--raw-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["raw"],
                    "help": "File suffix for raw input sequences (default: %(default)s).",
                },
            },
            {
                "args": ["--decompress-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["decompress"],
                    "help": "File suffix for decompressed sequences (default: %(default)s).",
                },
            },
            {
                "args": ["--merge-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["merge"],
                    "help": "File suffix for merged sequences (default: %(default)s).",
                },
            },
            {
                "args": ["--cutprimer-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["cutprimer"],
                    "help": "File suffix for trimmed sequences (default: %(default)s).",
                },
            },
            {
                "args": ["--dereplicate-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["dereplicate"],
                    "help": "File suffix for dereplicated sequences (default: %(default)s).",
                },
            },
            {
                "args": ["--denoise-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["denoise"],
                    "help": "File suffix for denoised sequences (default: %(default)s).",
                },
            },
            {
                "args": ["--assigntaxa-suffix"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["suffix"]["assigntaxa"],
                    "help": "File suffix for taxonomic assignment results (default: %(default)s).",
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
                "args": ["--min-read-len"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cutprimer"]["min_read_len"],
                    "help": "Minimum length of processed reads (default: %(default)s).",
                },
            },
            {
                "args": ["--max-read-len"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["cutprimer"]["max_read_len"],
                    "help": "Maximum length of processed reads (default: %(default)s).",
                },
            },
        ],
        "Denoising Settings": [
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
        "Taxonomic Assignment Settings": [
            {
                "args": ["--evalue"],
                "kwargs": {
                    "type": float,
                    "default": SETTINGS["assigntaxa"]["evalue"],
                    "help": "Expectation value (E) threshold for saving hits (default: %(default)s).",
                },
            },
            {
                "args": ["--qcov-hsp-perc"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["assigntaxa"]["qcov_hsp_perc"],
                    "help": "Percentage of query coverage for high-scoring segment pairs (default: %(default)s).",
                },
            },
            {
                "args": ["--perc-identity"],
                "kwargs": {
                    "type": int,
                    "default": SETTINGS["assigntaxa"]["perc_identity"],
                    "help": "Minimum Percentage identity for taxonomic assignment (default: %(default)s).",
                },
            },
            {
                "args": ["--specifiers"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["assigntaxa"]["specifiers"],
                    "help": "Custom format specifiers for BLAST results. (default: %(default)s).",
                },
            },
            {
                "args": ["-bdb", "--blast-db"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["assigntaxa"]["blast_db"],
                    "help": "Path to the BLAST database (default [NCBI remote nucleotide sequence database]: %(default)s).",
                },
            },
            {
                "args": ["-lndb", "--lineage-db"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["assigntaxa"]["lineage_db"],
                    "help": "Path to the lineage database (default [NCBI remote nucleotide lineage database]: %(default)s).",
                },
            },
            {
                "args": ["-email", "--entrez-email"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["assigntaxa"]["entrez_email"],
                    "help": "The email used by NCBI to contact you in case of excessive usage or issues. (default : %(default)s).",
                },
            },
        ],
        "External Program Settings": [
            {
                "args": ["--usearch-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["usearch"],
                    "help": "Command to execute USEARCH for merge, dereplicate, and denoise stages (default: %(default)s).",
                },
            },
            {
                "args": ["--cutadapt-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["cutadapt"],
                    "help": "Command to execute Cutadapt for primer trimming stage (default: %(default)s).",
                },
            },
            {
                "args": ["--blast-prog"],
                "kwargs": {
                    "type": str,
                    "default": SETTINGS["prog"]["blast"],
                    "help": "Command to execute BLAST for taxonomic assignment stage (default: %(default)s).",
                },
            },
        ],
        "Configuration Settings": [
            {
                "args": ["--verbose"],
                "kwargs": {
                    "action": "store_true",
                    "default": SETTINGS["config_basic"]["verbose"],
                    "help": "Enable detailed logging output (default: %(default)s).",
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

    options = vars(parser.parse_args())

    from ..bp import BioPipeline

    BioPipeline(**options)


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
