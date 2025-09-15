from ..common import base_logger

SETTINGS = {
    "enabled_stages": [
        "decompress",
        "merge",
        "cutprimer",
        "fqtofa",
        "dereplicate",
        "denoise",
        "assigntaxa",
    ],
    "dir_name": {
        "decompress": "decompress",
        "merge": "merge",
        "cutprimer": "cutprimer",
        "fqtofa": "fqtofa",
        "dereplicate": "dereplicate",
        "denoise": "denoise",
        "assigntaxa": "assigntaxa",
    },
    "suffix": {
        "raw": "_R1.fastq.gz",
        "decompress": "_R1.fastq",
        "merge": "_merged.fastq",
        "cutprimer": "_trimmed.fastq",
        "dereplicate": "_uniqs.fasta",
        "denoise": "_zotus.fasta",
        "assigntaxa": "_taxa.csv",
        "report": "_report.txt",
    },
    "merge": {
        "maxdiff": 5,
        "pctid": 90,
    },
    "cutprimer": {
        "rm_p_5": "GTCGGTAAAACTCGTGCCAGC",
        "rm_p_3": "CAAACTGGGATTAGATACCCCACTATG",
        "error_rate": 0.15,
        "min_read_len": 163,
        "max_read_len": 185,
    },
    "denoise": {
        "minsize": 8,
        "alpha": 2,
    },
    "assigntaxa": {
        "evalue": 0.00001,
        "qcov_hsp_perc": 90,
        "perc_identity": 90,
        "specifiers": "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "blast_db": "nt",
        "lineage_db": "nucleotide",
        "entrez_email": None,
    },
    "prog": {
        "usearch": "usearch",
        "cutadapt": "cutadapt",
        "blast": "blastn",
    },
    "config_basic": {
        "logger": base_logger.logger,
        "verbose": False,
        "dry": False,
    },
    "config_machine": {
        "n_cpu": 1,
        "memory": 8,
    },
}
