# eDNA Bioinformatics Pipeline (ednabp)

A pipeline for processing environmental DNA (eDNA) sequences through bioinformatics workflows including quality control, taxonomic assignment, and diversity analysis.

## Contents

- [Installation](#installation)
- [Modules](#modules)
- [Usage Examples](#usage-examples)
- [Testing](#testing)

## Installation


### Package Installation

1. Download from source:
```bash
git clone https://github.com/ComputationalAgronomy/ednabp.git
```
```bash
cd ednabp
````

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Install package locally:
```bash
pip install -e .
```

### Additional prerequisites for running the `bp` module

If you will run the `bp` module, ensure the following external tools are installed and available in your PATH:
* **Cutadapt** - [Download](https://cutadapt.readthedocs.io/en/stable/installation.html)
* **USEARCH** - [Download](https://github.com/rcedgar/usearch12.git)
* **NCBI-BLAST+** - [Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Additional prerequisites for running the `seq` module

If you will run the `seq` module, ensure the following external tools are installed and available in your PATH:
* **USEARCH** - [Download](https://github.com/rcedgar/usearch12.git)
* **Clustal Omega** - [Download](http://www.clustal.org/omega/)
* **IQTREE** - [Download](https://iqtree.github.io/doc/Quickstart#installation)

## Modules

![modules](./img/modules.png)

### 1. Bioinformatics Pipeline (`bp`)
Core processing pipeline with the following stages:

- **Decompress**: Extract compressed FASTQ files
- **Merge**: Combine paired-end reads
- **Cut Primer**: Remove primer sequences and length filtering
- **FASTQ to FASTA**: Format conversion
- **Dereplicate**: Remove duplicate sequences
- **Denoise**: Error correction
- **BLAST**: Taxonomic classification against a reference database
- **Add Lineage**: Annotate BLAST hits with full taxonomic lineage
- **LCA**: Resolve multiple hits per query to a consensus taxonomy using Lowest Common Ancestor
- **Add Haplotype**: Enrich results with haplotype sequence and abundance data

### 2. Data Management (`data`)

- **Data Objects**: Structured data containers for pipeline results.

A complete data container structure looks like the following:

![data](./img/data.png)

### 3. Diversity Analysis (`div`)

- **Writing**: Export diversity metrics CSV tables.
- **Plotting**: Visualization tools (`barchart`, `heatmap`, `rankcorr`, `sankey`) using Plotly as the underlying package.

### 4. Sequence Analysis (`seq`)

- **Clustering**: Sequence clustering analysis architecture that accepts a reducer class (e.g., `PCA`, `TSNE`, `UMAP`) and a clusterer class (e.g., `AgglomerativeClustering`, `HDBSCAN`). Note: You may need to install additional packages to access these classes.
- **Phylogenetics**: Tree construction and analysis using **IQTREE**.        .
- (TODO) **Haplotype Networks**: Write NEXUS files as input for **POPART** to draw haplotype networks.

May separate the `cluster` module into an independent repository in the future to keep each repo simple.

May remove the `phylo` and `hap_net modules` as it is somewhat redundant to use a Python interface rather than using those software packages directly.

## Usage Examples
- [Pipeline processing](#pipeline-processing)
- [Data management](#data-management)
- [Diversity metrics summary](#diversity-metrics-summary)

### Pipeline Processing

The BioPipeline processes sequencing data through a series of quality control and analysis stages. It expects input files in a specific format and produces organized output with intermediate results at each processing stage.

#### Input Format

- **Input directory**: Contains compressed FASTQ files (`.fastq.gz`)
- **For paired-end reads**: Files must be named with `_R1` and `_R2` suffixes (e.g., `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`)
- **For reads not requiring merge**: Files without read markers (e.g., `sample1.fastq.gz`); applicable to single-end, long-read, or already-merged data

#### Output Structure

```
output_path/
├── decompress/              # Decompressed FASTQ files
│   ├── sample1_R1.fastq
│   ├── sample1_R2.fastq
│   └── ...
├── merge/                   # Merged paired-end reads (or copies for single-end)
│   ├── sample1.fastq
│   └── ...
├── cutprimer/               # Primer-trimmed and length-filtered sequences
│   ├── sample1.fastq
│   └── ...
├── fqtofa/                  # Converted to FASTA format
│   ├── sample1.fasta
│   └── ...
├── dereplicate/             # Deduplicated sequences
│   ├── sample1.fasta
│   └── ...
├── denoise/                 # Denoised sequences (ZOTUs)
│   ├── sample1.fasta
│   ├── sample1_denoise_report.txt
│   └── ...
├── blast/                   # Taxonomic assignments
│   ├── sample1.csv          # Sequentially enriched: BLAST → lineage → LCA → haplotype
│   └── ...
└── stages.log               # Processing log
```

#### Run default pipeline

```python
from ednabp.bp import BioPipeline

# Process paired-end reads from directory containing _R1.fastq.gz and _R2.fastq.gz files
pipeline = BioPipeline(
    input_path="/path/to/files_folder",   # Directory containing multiple .fastq.gz files
    # input_path="/path/to/single_file",  # Alternative: single file input
    output_path="/path/to/output",
)
```

The pipeline automatically:
1. Decompresses `.fastq.gz` files
2. Merges paired-end reads (R1 + R2)
3. Removes primer sequences and filters by read length
4. Converts FASTQ to FASTA format
5. Removes duplicate sequences
6. Performs error correction (denoising)
7. Assigns taxonomy using BLAST (top 20 hits per query by default)
8. Adds full taxonomic lineage to each BLAST hit
9. Resolves multiple hits per query to a consensus taxonomy via LCA
10. Enriches with haplotype sequence and abundance data

#### Run custom settings

```python
custom_settings = {
    # Primer sequences (customize for your experiment)
    "rm_p_5": "GGACGATAAGACCCTATAAA",
    "rm_p_3": "ACTTTAGGGATAACAGCGT",
    
    # Length filtering parameters
    "min_read_len": 154,
    "max_read_len": 189,
    
    # Blast settings
    "blast_db": "/path/to/custom/blast/db",
    "lineage_db": "/path/to/custom/lineage/db",
    "maxhitnum": 100,
    
    # LCA settings
    "tol_pct": 1.0,
    
    # Performance settings
    "verbose": True,
    "n_cpu": 8,
}

pipeline = BioPipeline(
    input_path="/path/to/files_folder",
    output_path="/path/to/output",
    **custom_settings
)
```

#### Run partial pipeline

You can start the pipeline at any stage by specifying `enabled_stages`:

```python
# Start from merge stage (skip decompress)
settings = {
    "enabled_stages": ["merge", "cutprimer", "fqtofa", "dereplicate", "denoise", "blast", "addlineage", "lca", "addhap"],
}
pipeline = BioPipeline(
    input_path="/path/to/files_with_R1_R2_fastq",  # decompressed FASTQ files
    output_path="/path/to/output",
    **settings
)

# Start from fqtofa stage (skip decompress and merge)
settings = {
    "enabled_stages": ["fqtofa", "dereplicate", "denoise", "blast", "addlineage", "lca", "addhap"],
}
pipeline = BioPipeline(
    input_path="/path/to/files_with_merged_trimmed_fastq",
    output_path="/path/to/output",
    **settings
)
```
#### CLI

```bash
ednabp -i INPUT_PATH -o OUTPUT_PATH
```

To check the parameters, please run command:
```bash
ednabp -h
```

### Data Management
#### Import from ednabp.bp.BioPipeline outputs
```python
from ednabp.data import BPData

data = BPData()
data.import_data("results/")
# optional
data.import_metadata("path/to/sample_metadata")
data.import_spc_info("path/to/fishbase_db", "path/to/stock_db")
```

#### Import from MiFish Pipeline outputs

This package supports import data from another popular pipeline to run downstream analysis.

[MiFish Pipeline webpage](https://mitofish.aori.u-tokyo.ac.jp/mifish/)

```python
from ednabp.data import MitoData

data = MitoData()
data.import_data("results/")
# optional
data.import_metadata("path/to/sample_metadata")
```

#### Reuse a data container
You can serialize and deserialize a data container for repetitive use. This process is known as "pickling" and "unpickling." Note: Never unpickle a .pkl file from an unknown source.

```python
data.pickle_data("path/to/save_dir", "save_name")
```
Next time, you only need to unpickle the data container and don't need to import everything again.

```python
data = BPData() # or MitoData()
data.unpickle_data("path/to/pkl_file")
```

### Diversity Metrics Summary
Here is an example with writing abundance table and drawing barchart of species abundance across samples.


```python
from ednabp.div.write import Writer
from ednabp.div.plot import barchart

# Create abundance dataframe
writer = Writer(data)
df = writer.abundance(taxa_lv='species')

# Generate barchart
fig, plotter = barchart(
    df=df,
    values='abundance',
    index='species',
    columns='sample_id'
)
```

We also provide two other metrics: **richness** and **detection probability**, plus three additional visualization options: **heatmap**, **sankey diagram**, and **rank correlation matrix**. Additionally, you can customize parameters for summarizing metrics and visualizing data, such as `taxa_lv`, `values`, `index`, and `columns`. These options give you the flexibility to describe your own data.

## Testing

Run the test suite:

```bash
pytest
# or
pytest ./tests/test_XXX.py
```
