# eDNA_Bioinformatics_Pipeline
The purpose of this project is to implement eDNA bioinformatics processing and to explore taxonomic delimitation from various perspectives using eDNA sequences.

## Contents

- [Getting Started](#getting-started)
  - [Docker Version](#docker-version)
    - [Useful Docker Commands](#useful-docker-commands)
  - [Local Version](#local-version)
- [Pytest](#pytest)
- [Usage](#usage)
  - [`fastq_processor`](#fastq_processor-module)
    - [Simplist Example](#simplist-example)
    - [Optional Parameters](#other-optional-parameters)
  - [`analysis_toolkit`](#analysis_toolkit-module)

## Getting Started

### Docker Version

**Step 1.** Download [Docker Desktop](https://docs.docker.com/get-docker/) on your own machine.

After installation, run `docker version` in terminal. If the version is displayed, it means the installation was successful.

**Step 2.** Clone the GitHub repo:
```sh
git clone https://github.com/ComputationalAgronomy/eDNA_bioinformatics_pipeline.git
```

**Step 3.** Build an image according to the Dockerfile.
```sh
docker build -t [ImageName] .
```

After the Dockerfile is successfully exported to an image, the [ImageName] should appear in the repo list if you use `docker image ls` to check.

**Step 4.** launch a new container from the Docker image that was just built:
```sh
docker run -it [ImageName]
```

If the launch is successful, your terminal should display something like `(base) root@93f4d3cf355f:/#`.

`(base)` indicates the conda environment you are using, where all dependent Python packages are installed (don't deactivate this!). `93f4d3cf355f` indicates the ID of this container.

**Step fin.** The Container is Ready to Work! Let's try [the first example](#simplist-example)!

#### Useful Docker Commands

(base) root@93f4d3cf355f:/# `exit` or `Ctrl+Z`: Exit the container.

`docker cp [container ID]:/path/to/file /host/destination/folder`: Copying files from Docker container to host.

`docker cp /path/to/file [container ID]:/container/destination/folder`: Copying files from host to Docker container.

`docker exec -it [container ID or Name] bash`: Enter a running container's shell.

`docker container ls -a`: List all containers

`docker stop [container ID or Name]`: Stop a running container.

`docker start [container ID or Name]`: Start a stopped container.

`docker rmi [ImageName]`: Remove a Docker image.

`docker container rm [container ID or Name]`: Remove a container.

`docker system prune (--force)`: Remove \<none> TAG images (be careful when using this command).

### Local Version

#### Dependency Installation

Make sure you have installed all of the following prerequisites on your machine:
* BBMap - [Download](https://sourceforge.net/projects/bbmap/) (You may also need to install java if you haven't already)
* Clustal Omega - [Download](http://www.clustal.org/omega/)
* Cutadapt - [Download](https://cutadapt.readthedocs.io/en/stable/installation.html)
* IQTREE2 - [Download](http://www.iqtree.org/)
* NCBI-BLAST+ - [Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* USEARCH - [Download](https://www.drive5.com/usearch/download.html)

, and ensure the path of downloaded software is added to the "Path" variable in "Environmental Variables".

#### Installation
1. Required Python Package Installation
```sh
python -m pip install -r requirements.txt
```
2. Local package installs
```sh
python -m pip install -e .
```

## Pytest
check `pytest.ini`
```sh
pytest
# OR
pytest tests
```

## Usage

### `fastq_processor` Module
This module processes raw FASTQ data through several stages including paired-end merging, primer cutting, reformatting, dereplication, denoising, and taxonomic assignment.

#### Simplist Example
**Step 0.** By default, your container will start in the `/workplace` folder. Use the `tree` command to check the prepared materials for this practice!

**Step 1.** Use a text editor to create a Python script.
```sh
nano [FileName].py
```

**Step 2.** Write the following code, then save and exit the text editor.
```python
from fastq_processor import FastqProcessor

FastqProcessor(
    stages_parent_dir="./fastq_processor/stages",
    fastq_dir_name="fastq",
    db_path="./fastq_processor/db/MiFish",
    lineage_path="./fastq_processor/lineage/lineage.csv",
)
```

**Ensure that `stages_parent_dir` exists and contains a subdirectory named `fastq_dir_name` with your raw FASTQ files.**
```
stages/
└── fastq/
    ├── sample1_R1.fastq
    └── sample1_R2.fastq
```

The **`db_path`** serves as the reference data for assigning denoised sequences to species. **It should be set to the folder path containing the indexed files, with the prefix string added.** The MiFish index files provided in the *example* folder are built using the [complete+partial mtDNA sequence file](https://mitofish.aori.u-tokyo.ac.jp/species/detail/download/?filename=download%2F/complete_partial_mitogenomes.zip) downloaded from the MiFish Pipeline. If you want to use a custom FASTA file as the reference database, you need to create the index using the `makeblastdb` command from ncbi-blast+. Here is the command: `makeblastdb -in ref.fasta -dbtype nucl -out db_prefix`

The **`lineage_path`** serves as a reference for labeled species at the taxonomic level above the species. **It should be set to the path of the lineage file.** The [lineage.csv](https://github.com/billzt/MiFish/blob/main/mifish/data/lineage.csv) provided in the *example* folder is downloaded from the MiFish GitHub repository. If you want to use another lineage file, please ensure it includes taxonomic information from the domain to genus level and maintains the same format.

**Step 3.** Run the Python script you just wrote in your container.
```sh
python [FileName].py
```

**Step 4.** Use the `tree` command again to check the results. It should look something like this:
```
stages/
├── fastq/
|   ├── sample1_R1.fastq
|   └── sample1_R2.fastq
├── merge/
|   ├── sample1_merge.fastq
|   └── sample1_report.txt
├── cut_primer/
|   ├── sample1_cut.fastq
|   └── sample1_report.txt
├── fq_to_fa/
|   ├── sample1_cut.fasta
|   └── sample1_report.txt
├── dereplicate/
|   ├── sample1_uniq.fasta
|   └── sample1_report.txt
├── denoise/
|   ├── sample1_denoise.fasta
|   ├── sample1_denoise_report.txt
|   └── sample1_report.txt
└── denoise/
    └── sample1_blast.csv
```

**Step 5.** Copy the result files from the Docker container back to the host. 

**In the host terminal**, use the following command:
```sh
docker cp [container ID]:/path/to/stages/folder /host/destination/folder
```
**You have successfully run the process and obtained the processed data!**

#### Optional Parameters:

`enabled_stages`: The list of stages to run.
 Default is `["merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "assigntaxa"]`.

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

### `Analysis_toolkit` Module