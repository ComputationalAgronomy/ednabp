#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from .fastq_processor.run_processor import FastqProcessor

def main():
    parser = argparse.ArgumentParser(
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS,
        description="A pipeline for processing eDNA bioinformatics workflows."
    )
    
    # Required Input & Output Path
    io_group = parser.add_argument_group("Required Input & Output Path")
    io_group.add_argument("-i", "--input-path", required=True, help="Input FASTA sequence file path or directory.")
    io_group.add_argument("-o", "--output-path", required=True, help="Directory where output files will be saved.")
    
    # Stages Executed
    stage_group = parser.add_argument_group("Stages Executed")
    stage_group.add_argument("-e", "--enabled-stages", nargs="+", help="Process stages to be executed in order (default [run all stages]: 'decompress merge cutprimer fqtofa dereplicate denoise assigntaxa').")
    
    # Directory Names
    dir_group = parser.add_argument_group("Directory Names")
    dir_group.add_argument("--decompress-dir-name", help="Name of the subdirectory for the decompression stage (default: 'decompress').")
    dir_group.add_argument("--merge-dir-name", help="Name of the subdirectory for the merge stage (default: 'merge').")
    dir_group.add_argument("--cutprimer-dir-name", help="Name of the subdirectory for the cut-primer stage (default: 'cutprimer').")
    dir_group.add_argument("--fqtofa-dir-name", help="Name of the subdirectory for FastQ to FastA conversion (default: 'fqtofa').")
    dir_group.add_argument("--dereplicate-dir-name", help="Name of the subdirectory for the dereplication stage (default: 'dereplicate').")
    dir_group.add_argument("--denoise-dir-name", help="Name of the subdirectory for the denoising stage (default: 'denoise').")
    dir_group.add_argument("--assigntaxa-dir-name", help="Name of the subdirectory for taxonomic assignment (default: 'assigntaxa').")
    
    # File Suffixes
    suffix_group = parser.add_argument_group("File Suffixes")
    suffix_group.add_argument("--raw-suffix", help="File suffix for raw input sequences (default: '_R1.fastq.gz').")
    suffix_group.add_argument("--decompress-suffix", help="File suffix for decompressed sequences (default: '_R1.fastq').")
    suffix_group.add_argument("--merge-suffix", help="File suffix for merged sequences (default: '_merged.fastq').")
    suffix_group.add_argument("--cutprimer-suffix", help="File suffix for trimmed sequences (default: '_trimmed.fastq').")
    suffix_group.add_argument("--dereplicate-suffix", help="File suffix for dereplicated sequences (default: '_uniqs.fasta').")
    suffix_group.add_argument("--denoise-suffix", help="File suffix for denoised sequences (default: '_zotus.fasta').")
    suffix_group.add_argument("--assigntaxa-suffix", help="File suffix for taxonomic assignment results (default: '_taxa.csv').")
    
    # Merge Settings
    merge_group = parser.add_argument_group("Paired-end Merge Settings")
    merge_group.add_argument("--maxdiff", type=int, help="Maximum number of mismatches in alignment during merging (default: 5).")
    merge_group.add_argument("--pctid", type=int, help="Minimum percentage identity required for merging reads (default: 90).")
    
    # Primer Cutting Settings
    cuter_group = parser.add_argument_group("Primer Trimming Settings")
    cuter_group.add_argument("--rm-p-5", help="Non-internal 5’ primer (default: 'GTCGGTAAAACTCGTGCCAGC').")
    cuter_group.add_argument("--rm-p-3", help="Non-internal 3’ primer (default: 'CAAACTGGGATTAGATACCCCACTATG').")
    cuter_group.add_argument("--error-rate", type=float, help="Maximum allowable error rate for primer matching (default: 0.15).")
    cuter_group.add_argument("--min-read-len", type=int, help="Minimum length of processed reads (default: 204).")
    cuter_group.add_argument("--max-read-len", type=int, help="Maximum length of processed reads (default: 254).")
    
    # Denoising Settings
    denoise_group = parser.add_argument_group("Denoising Settings")
    denoise_group.add_argument("--minsize", type=int, help="Minimum abundance of sequences to retain during denoising (default: 8).")
    denoise_group.add_argument("--alpha", type=int, help="Denoising sensitivity parameter (default: 2).")
    
    # Taxonomic Assignment Settings
    taxa_group = parser.add_argument_group("Taxonomic Assignment Settings")
    taxa_group.add_argument("-db", "--db-path", help="Path to the taxonomic database (default: None).")
    taxa_group.add_argument("-lng", "--lineage-path", help="Path to the taxonomic lineage file (default: None).")
    taxa_group.add_argument("--evalue", type=float, help="Expectation value (E) threshold for saving hits (default: 0.00001).")
    taxa_group.add_argument("--qcov-hsp-perc", type=int, help="Percentage of query coverage for high-scoring segment pairs (default: 90).")
    taxa_group.add_argument("--perc-identity", type=int, help="Minimum percentage identity for taxonomic assignment (default: 90).")
    taxa_group.add_argument("--specifiers", help="Custom format specifiers for BLAST results. (default: 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore').")
    
    # Configuration Settings
    config_group = parser.add_argument_group("Configuration Settings")
    config_group.add_argument("--verbose", action="store_true", help="Enable detailed logging output (default: True).")
    config_group.add_argument("--n-cpu", type=int, help="Number of CPU cores to use for processing (default: 1).")

    options = vars(parser.parse_args())

    FastqProcessor(**options)

if __name__ == "__main__":
    main()
