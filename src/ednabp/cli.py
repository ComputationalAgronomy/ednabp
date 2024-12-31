#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from .fastq_processor.run_processor import FastqProcessor


def main():
    argument_groups = {
        "Required Input & Output Path": [
            {"args": ["-i", "--input-path"], "kwargs": {"required": True, "help": "Input FASTA sequence file path or folder path."}},
            {"args": ["-o", "--output-path"], "kwargs": {"required": True, "help": "Directory where output files will be saved."}}
        ],
        "Stages Executed": [
            {"args": ["-e", "--enabled-stages"], "kwargs": {"nargs": "+", "help": "Process stages to be executed in order (default [run all stages]: decompress merge cutprimer fqtofa dereplicate denoise assigntaxa)."}}
        ],
        "Directory Names": [
            {"args": ["--decompress-dir-name"], "kwargs": {"help": "Name of the subdirectory for the decompression stage (default: decompress)."}},
            {"args": ["--merge-dir-name"], "kwargs": {"help": "Name of the subdirectory for the merge stage (default: merge)."}},
            {"args": ["--cutprimer-dir-name"], "kwargs": {"help": "Name of the subdirectory for the cut-primer stage (default: cutprimer)."}},
            {"args": ["--fqtofa-dir-name"], "kwargs": {"help": "Name of the subdirectory for FastQ to FastA conversion (default: fqtofa)."}},
            {"args": ["--dereplicate-dir-name"], "kwargs": {"help": "Name of the subdirectory for the dereplication stage (default: dereplicate)."}},
            {"args": ["--denoise-dir-name"], "kwargs": {"help": "Name of the subdirectory for the denoising stage (default: denoise)."}},
            {"args": ["--assigntaxa-dir-name"], "kwargs": {"help": "Name of the subdirectory for taxonomic assignment (default: assigntaxa)."}}
        ],
        "File Suffixes": [
            {"args": ["--raw-suffix"], "kwargs": {"help": "File suffix for raw input sequences (default: _R1.fastq.gz)."}},
            {"args": ["--decompress-suffix"], "kwargs": {"help": "File suffix for decompressed sequences (default: _R1.fastq)."}},
            {"args": ["--merge-suffix"], "kwargs": {"help": "File suffix for merged sequences (default: _merged.fastq)."}},
            {"args": ["--cutprimer-suffix"], "kwargs": {"help": "File suffix for trimmed sequences (default: _trimmed.fastq)."}},
            {"args": ["--dereplicate-suffix"], "kwargs": {"help": "File suffix for dereplicated sequences (default: _uniqs.fasta)."}},
            {"args": ["--denoise-suffix"], "kwargs": {"help": "File suffix for denoised sequences (default: _zotus.fasta)."}},
            {"args": ["--assigntaxa-suffix"], "kwargs": {"help": "File suffix for taxonomic assignment results (default: _taxa.csv)."}}
        ],
        "Paired-end Merge Settings": [
            {"args": ["--maxdiff"], "kwargs": {"type": int, "help": "Maximum number of mismatches in alignment during merging (default: 5)."}},
            {"args": ["--pctid"], "kwargs": {"type": int, "help": "Minimum percentage identity required for merging reads (default: 90)."}}
        ],
        "Primer Trimming Settings": [
            {"args": ["--rm-p-5"], "kwargs": {"help": "Non-internal 5’ primer (default [MiFish-UF]: GTCGGTAAAACTCGTGCCAGC)."}},
            {"args": ["--rm-p-3"], "kwargs": {"help": "Non-internal 3’ primer (default [rev-com MiFish-UR]: CAAACTGGGATTAGATACCCCACTATG)."}},
            {"args": ["--error-rate"], "kwargs": {"type": float, "help": "Maximum allowable error rate for primer matching (default: 0.15)."}},
            {"args": ["--min-read-len"], "kwargs": {"type": int, "help": "Minimum length of processed reads (default: 204)."}},
            {"args": ["--max-read-len"], "kwargs": {"type": int, "help": "Maximum length of processed reads (default: 254)."}}
        ],
        "Denoising Settings": [
            {"args": ["--minsize"], "kwargs": {"type": int, "help": "Minimum abundance of sequences to retain during denoising (default: 8)."}},
            {"args": ["--alpha"], "kwargs": {"type": int, "help": "Denoising sensitivity parameter (default: 2)."}}
        ],
        "Taxonomic Assignment Settings": [
            {"args": ["-db", "--db-path"], "kwargs": {"help": "Path to the taxonomic database (default: None)."}},
            {"args": ["-ln", "--lineage-path"], "kwargs": {"help": "Path to the taxonomic lineage file (default: None)."}},
            {"args": ["--evalue"], "kwargs": {"type": float, "help": "Expectation value (E) threshold for saving hits (default: 0.00001)."}},
            {"args": ["--qcov-hsp-perc"], "kwargs": {"type": int, "help": "Percentage of query coverage for high-scoring segment pairs (default: 90)."}},
            {"args": ["--perc-identity"], "kwargs": {"type": int, "help": "Minimum Percentage identity for taxonomic assignment (default: 90)."}},
            {"args": ["--specifiers"], "kwargs": {"help": "Custom format specifiers for BLAST results. (default: 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore')."}}
        ],
        "Configuration Settings": [
            {"args": ["--verbose"], "kwargs": {"action": "store_true", "help": "Enable detailed logging output (default: True)."}},
            {"args": ["--n-cpu"], "kwargs": {"type": int, "help": "Number of CPU cores to use for processing (default: 1)."}}
        ]
    }

    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: MyFormatter(prog, max_help_position=32),
        argument_default=argparse.SUPPRESS,
        description = "A pipeline for processing eDNA bioinformatics workflows."
    )

    for group_name, args_list in argument_groups.items():
        group = parser.add_argument_group(group_name)
        for arg_data in args_list:
            group.add_argument(*arg_data["args"], **arg_data["kwargs"])

    options = vars(parser.parse_args())

    FastqProcessor(**options)

class MyFormatter(argparse.HelpFormatter):
    """
        for matt wilkie on SO
        https://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """
    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
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
                args_string = self._format_args(action, '')
                for option_string in action.option_strings:
                    parts.append(option_string)
                return '%s %s' % (', '.join(parts), args_string)
            return ', '.join(parts)

if __name__ == "__main__":
    main()
