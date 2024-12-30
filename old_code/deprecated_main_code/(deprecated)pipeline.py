import os
import pandas as pd
import csv
import numpy as np
import subprocess
# import duckdb


def get_prefix_with_suffix(in_dir, suffix):
    files = os.listdir(in_dir)
    prefix = []
    for filename in files:
        if filename.endswith(suffix):
            prefix.append(filename.replace(suffix, ''))
    return prefix

def merge_fq(in_dir, out_dir, maxdiff=5, pctid=90, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_R1.fastq')
    for num in prefix:
        cmd = f'usearch -fastq_mergepairs {in_dir}{num}_R1.fastq -fastqout {out_dir}{num}_merge.fastq -threads {cpu} \
                -fastq_maxdiffs {maxdiff} -fastq_pctid {pctid} -report {out_dir}{num}_report.txt'
        subprocess.run(cmd, shell=True)

def cut_adapt(in_dir, out_dir, rm_p_5='GTCGGTAAAACTCGTGCCAGC', rm_p_3='CAAACTGGGATTAGATACCCCACTATG', min_read_len=204, max_read_len=254, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_merge.fastq')
    for num in prefix:
        cmd = (f'cutadapt -g "{rm_p_5};max_error_rate=0.15...{rm_p_3};max_error_rate=0.15" -j {cpu} \
                {in_dir}{num}_merge.fastq --discard-untrimmed \
                -m {min_read_len-len(rm_p_5)-len(rm_p_3)} -M {max_read_len-len(rm_p_5)-len(rm_p_3)} \
                >{out_dir}{num}_cut.fastq 2>{out_dir}{num}_report.txt')
        subprocess.run(cmd, shell=True)

def fq_to_fa(in_dir, out_dir, bbmap_dir):
    prefix = get_prefix_with_suffix(in_dir, '_cut.fastq')
    for num in prefix:
        cmd = f'bash {bbmap_dir}reformat.sh in={in_dir}{num}_cut.fastq out={out_dir}{num}_processed.fasta'
        subprocess.run(cmd, shell=True)

def dereplicate(in_dir, out_dir, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_processed.fasta')
    for num in prefix:
        cmd = f'usearch -fastx_uniques {in_dir}{num}_processed.fasta -threads {cpu} \
                -sizeout -relabel Uniq -fastaout {out_dir}{num}_derep.fasta \
                >{out_dir}{num}_report.txt 2>&1'
        subprocess.run(cmd, shell=True)

def cluster_otu(in_dir, out_dir, minsize=2, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_derep.fasta')
    for num in prefix:
        cmd = f'usearch -cluster_otus {in_dir}{num}_derep.fasta -minsize {minsize} -threads {cpu} \
                -otus {out_dir}{num}_otu.fasta -uparseout {out_dir}{num}_otu_report.txt -relabel Otu \
                >{out_dir}{num}_report.txt 2>&1'
        subprocess.run(cmd, shell=True)

def cluster_zotu(in_dir, out_dir, minsize=8, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_derep.fasta')
    for num in prefix:
        cmd = f'usearch -unoise3 {in_dir}{num}_derep.fasta -minsize {minsize} -threads {cpu} \
                -zotus {out_dir}{num}_zotu.fasta -tabbedout {out_dir}{num}_zotu_report.txt \
                >{out_dir}{num}_report.txt 2>&1'
        subprocess.run(cmd, shell=True)

def blast_otu(in_dir, out_dir, db_path, lineage_path, otu_type='otu', cpu=0, maxhitnum=5, specifiers='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'):
    if otu_type not in ['otu', 'zotu']:
        print("Error: 'otu_type' must be 'otu' or 'zotu'.")
        return

    prefix = get_prefix_with_suffix(in_dir, f'_{otu_type}.fasta')
    for num in prefix:
        cmd = f'blastn -db {db_path} -query {in_dir}{num}_{otu_type}.fasta -outfmt "10 {specifiers}" \
                -max_target_seqs {maxhitnum} -evalue 0.00001 -qcov_hsp_perc 90 -perc_identity 90 \
                -out {out_dir}{num}_{otu_type}.csv'
        cmd = cmd + f' -num_threads {cpu}' if cpu!=0 else cmd
        subprocess.run(cmd, shell=True)

    genus2taxonomy = {}
    with open(lineage_path+'lineage.csv') as in_handle:
        csvReader = csv.DictReader(in_handle)
        for row in csvReader:
            genus2taxonomy[row['genus_name']] = [row['kingdom_name'], row['phylum_name'], row['class_name'], \
                row['order_name'], row['family_name']]
    for num in prefix:
        filename = f'{out_dir}{num}_{otu_type}'
        _add_taxonomy(filename=filename, genus2taxonomy=genus2taxonomy)

@staticmethod
def _add_taxonomy(filename, genus2taxonomy):
    blast_result = pd.read_csv(f'{filename}.csv', header=None)
    taxa_matrix = []
    for sseqid in blast_result[1]:
        sacc, species = sseqid.split('|')[1], sseqid.split('|')[-1]
        species_firstname = species.split('_')[0]
        if species_firstname in genus2taxonomy.keys():
            genus = species_firstname
        else:
            genus = [genus_level for genus_level, other_levals in genus2taxonomy.items() if species_firstname in other_levals][0]
        taxa_matrix.append([sacc, species, genus, genus2taxonomy[genus][4], genus2taxonomy[genus][3], genus2taxonomy[genus][2], genus2taxonomy[genus][1], genus2taxonomy[genus][0]])
    taxa_matrix = np.array(taxa_matrix)
    sacc_list, species_list, genus_list, family_list, order_list, class_list, phylum_list, kingdom_list = taxa_matrix[:, 0].tolist(), taxa_matrix[:, 1].tolist(), taxa_matrix[:, 2].tolist(), taxa_matrix[:, 3].tolist(), taxa_matrix[:, 4].tolist(), taxa_matrix[:, 5].tolist(), taxa_matrix[:, 6].tolist(), taxa_matrix[:, 7].tolist()
    blast_result[1] = sacc_list
    blast_result.insert(2, '2', species_list)
    blast_result.insert(3, '3', genus_list)
    blast_result.insert(4, '4', family_list)
    blast_result.insert(5, '5', order_list)
    blast_result.insert(6, '6', class_list)
    blast_result.insert(7, '7', phylum_list)
    blast_result.insert(8, '8', kingdom_list)
    blast_result.to_csv(f'{filename}.csv', index=False, header=None)

if __name__ == '__main__':

    in_dir = './cleandata/5_haploid/zotu/'
    out_dir = './cleandata/test/mifish_partial/'


    # merge_fq(in_dir=in_dir, out_dir=out_dir, cpu=6)

    # cut_adapt(in_dir=in_dir, out_dir=out_dir, cpu=6)

    # bbmap_dir = './bbmap/'
    # fq_to_fa(bbmap_dir=bbmap_dir, in_dir=in_dir, out_dir=out_dir)

    # dereplicate(in_dir=in_dir, out_dir=out_dir, cpu=16)

    # cluster_otu(in_dir=in_dir, out_dir=out_dir, minsize=8, cpu=16)
    # cluster_zotu(in_dir=in_dir, out_dir=out_dir, minsize=8, cpu=16)

    # mifish_path = './database/mifish_partial/mifish'
    # lineage_path = './database/'
    # blast_otu(in_dir=in_dir, out_dir=out_dir, db_path=mifish_path, lineage_path=lineage_path, maxhitnum=1, otu_type='zotu', cpu=16)

    # # ncbi_path = 'nt -remote'
    # blast_otu(in_dir=in_dir, out_dir=out_dir, db_path=ncbi_path, out_name='ncbi', otu_type='zotu')
    # read_blast_csv(in_dir=in_dir, out_dir=out_dir)

    # for expanding species related information in future
    # conn = duckdb.connect()
    # fishbase_file = './data/species.parquet'
    # stock_file = './data/stocks.parquet'
    # link_fishbase = conn.from_parquet(fishbase_file)
    # link_stock = conn.from_parquet(stock_file)