# eDNA_Bioinformatics_Pipeline
The purpose of this project is to implement eDNA bioinformatics processing and to explore taxonomic delimitation from various perspectives using eDNA sequences.

# Contents

- [Getting Started](#getting-started)
  
- [Usage](#usage)

- [Pytest](#pytest)

# Getting Started
 - [Local Version](#local-version)

 - [Google Colab Version](#google-colab-version)

 - [Docker Version](#docker-version)
## Local Version

### Dependency Installation

Make sure you have installed all of the following prerequisites on your machine:
* Clustal Omega - [Download](http://www.clustal.org/omega/)
* Cutadapt - [Download](https://cutadapt.readthedocs.io/en/stable/installation.html)
* IQTREE2 - [Download](http://www.iqtree.org/)
* NCBI-BLAST+ - [Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* USEARCH12 - [Download](https://github.com/rcedgar/usearch12.git)

, and ensure the path of downloaded software is added to the "Path" variable in "Environmental Variables".

### Installation
1. Required Python Package Installation
```sh
python -m pip install -r requirements.txt
```
2. Local package installs
```sh
python -m pip install -e .
```

## Google Colab Version

For details check [ednabp.ipynb](https://colab.research.google.com/drive/1HTjt7VIZwWkEjDWdPMAJ3U8vSWy2KlVd?usp=sharing)

## Docker Version

For details check `docker.md`

# Usage
  - [`fastq_processor`](#fastq_processor-module)

  - [`analysis_toolkit`](#analysis_toolkit-module)

## `fastq_processor` Module
This module processes raw FASTQ data through several stages including paired-end merging, primer cutting, reformatting, dereplication, denoising, and taxonomic assignment.

### Step 0. Prepare necessary files.
 - A directory to save the output files, and ensure that it contains a subdirectory with your raw FASTQ.GZ files. Filenames for raw data should follow the pattern `<prefix><suffix>`. You can change the suffix pattern by modifying the `raw_suffix` parameter (see [Optional Parameters](#optional-parameters)).
 ```
parent_dir/
└── raw_dir/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    └── ...
```

 - A directory contains index files, which serve as the reference data for assigning sequences to species. If you want to create a custom set of index files as a reference database, use the following command:
 ```
makeblastdb -in ref.fasta -dbtype nucl -out db_prefix`
```

 - A lineage CSV file includes taxonomic information from the domain to genus level, which serves as a reference for labeled species at the taxonomic level above the species. 

### Step 1. Write the following code to create a Python script and run the script.
```python
from fastq_processor import FastqProcessor

FastqProcessor(stages_parent_dir = "/path/to/parent_dir",
               raw_dir_name = "RAW_DIR_NAME",
               db_path = "/path/to/db_dir/DB_PREFIX",
               lineage_path = "path/to/lineage.csv")
```
The `stages_parent_dir` should be set to the directory path that contains a subdirectory with raw data to be processed and where output files will be saved.

The `raw_dir_name` should be set to the name of the subdirectory that contains raw data.

The `db_path` should be set to the directory path containing the indexed files, with the prefix string appended to the folder path.

The `lineage_path` should be set to the path of the lineage CVS file.

### Step 2. Check the results. It should look something like this:
```
parent_dir/
├── raw_dir/
|   ├── sample1_R1.fastq.gz
|   └── sample1_R2.fastq.gz
├── decompress/
|   ├── sample1_R1.fastq
|   └── sample1_R2.fastq
├── merge/
|   ├── sample1_merge.fastq
|   └── sample1_report.txt
├── cut_primer/
|   ├── sample1_cut.fastq
|   └── sample1_report.txt
├── fq_to_fa/
|   └── sample1_cut.fasta
├── dereplicate/
|   ├── sample1_uniq.fasta
|   └── sample1_report.txt
├── denoise/
|   ├── sample1_denoise.fasta
|   ├── sample1_denoise_report.txt
|   └── sample1_report.txt
└── blast/
    └── sample1_blast.csv
```

**You have successfully run the process and obtained the processed data!**

### Optional Parameters:

`enabled_stages`: The list of stages to run.
 Default is `["decompress", "merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "assigntaxa"]`.

`raw_suffix`: The suffix pattern used for R1 raw data. Default is `_R1.fastq.gz`.

`n_cpu`: Number of CPU cores to be used for processing. Default is `1`.

`verbose`: Set to True to enable detailed logging output. Default is `True`.

Paired-end merging related:

`maxdiff`: Maximum number of mismatches in the alignment. Default is `5`.

`pctid`: Minimum %id of alignment. Default is `90`.

Primer cutting related:

`rm_p_5`: Non-internal 5’ primer. Default is `"GTCGGTAAAACTCGTGCCAGC"` (MiFish-UF).

`rm_p_3`: Non-internal 3’ primer. Default is `"CAAACTGGGATTAGATACCCCACTATG"` (reverse-complement MiFish-UR).

`error_rate`: The maximum rate of error could be tolerated. The actual error rate is computed as the number of errors in the match divided by the length of the matching part of the primer. Default is `0.15`.

`min_read_len`: Discard processed reads that are shorter than this parameter. Default is `204`.

`max_read_len`: Discard processed reads that are longer than this parameter. Default is `254`.

Denoising related:

`minsize`: Discard sequences with abundance that are smaller than this parameter. Default is `8`.

`alpha`: See [UNOISE2 paper](https://www.biorxiv.org/content/10.1101/081257v1.full) for definition. Default is `2`.

Taxonomic assignment related:

`evalue`: Expectation value (E) threshold for saving hits. Default is `0.00001`.

`qcov_hsp_perc`: The percent threshold of the query sequence that has to form an alignment against the reference to be retained. Default is `90`.

`perc_identity`: Percent identity cutoff. Default is `90`.

`specifiers`: Use to customize format specifiers. Default is `"qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"`

## `Analysis_toolkit` Module
This module conducts various types of writing and plotting tasks using the output results of `fastq_processor`.

[Step 0. Create Data Object](#create-data-object)

Haplotypes related:

 - [UMAP]

 - [HDBSCAN on UMAP embedding]

 - [Maximum Likelihood tree]

 - [NEXUS file for Haplotype Network]

Abundance or Richness related:

 - [Barchart](#barchart)

 - [Heatmap](#heatmap)

### Create Data Object
#### Import Data
```python
from analysis_toolkit import SampleData

sample_data = SampleData()
sample_data.import_data(uniq_dir = "/path/to/dir_save_uniq_fa",
                        denoise_dir = "/path/to/dir_save_denoise_fa",
                        denoise_report_dir = "/path/to/dir_save_denoise_report",
                        blast_dir = "/path/to/dir_save_blast_csv",
                        sample_info_path = "/path/to/sample_info.csv")
```
#### Merge Data
```python
# Suppose you import another dataset into an object called "sample_data2",
# and you want to combine these two datasets to run an analysis.
sample_data.merge_data(sample_data2)
```
#### Save Instance
```python
sample_data.save_data(save_instance_dir = "/path/to/save",
                      save_prefix = "INSTANCE_NAME",
                      overwrite = True)
```

#### Load Data
```python
from analysis_toolkit import SampleData

sample_data = SampleData()
sample_data.load_data(load_instance_path = "/path/to/save/INSTANCE_NAME.pkl")
```

### Barchart
```python
from analysis_toolkit import BarchartRunner

br = BarchartRunner(sample_data)
br.run_write(taxa_level = "family",
             write_type = "abundance",
             save_dir = "stage_test",
             normalize = True)
br.run_plot(csv_path = "stage_test/species_abundance.csv")
```

### Heatmap
```python
from analysis_toolkit import HeatmapRunner

hr = HeatmapRunner(sample_data)
hr.run_write(taxa_level = "family",
             write_type = "richness",
             save_dir = "TI_test")
hr.run_plot(csv_path = "TI_test/Species_richness.csv",
            x_categories = ["Site", "Sample"])
```
## Pytest
check `pytest.ini`
```sh
pytest
# OR
pytest tests
```