import os
from tempfile import TemporaryDirectory
from typing import Literal

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.spatial.distance import pdist, squareform

from ..common import base_logger, config
from ..data import BPData, MitoData
from . import write


def encode_pdist_mx(
    fa_path,
):
    base_to_num = {"A": 1, "C": 2, "G": 3, "T": 4}
    seqs = []
    descriptions = []
    with open(fa_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            numeric_seq = [base_to_num.get(base, 0) for base in record.seq]
            seqs.append(numeric_seq)
            descriptions.append(record.description)

    seqs = np.array(seqs)
    distances = pdist(seqs, "hamming")
    matrix = squareform(distances)

    return matrix, descriptions


def encode_onehot_mx(fa_path):
    base_map = {
        "A": np.array([1, 0, 0, 0]),
        "C": np.array([0, 1, 0, 0]),
        "G": np.array([0, 0, 1, 0]),
        "T": np.array([0, 0, 0, 1]),
        "-": np.array([-1, -1, -1, -1]),
        "N": np.array([0.25, 0.25, 0.25, 0.25]),
    }

    sequences = []
    descriptions = []
    with open(fa_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
            descriptions.append(record.description)

    seq_len = len(sequences[0])
    matrix = np.zeros((len(sequences), seq_len * 4))

    for i, seq in enumerate(sequences):
        for j, base in enumerate(seq):
            matrix[i, j * 4 : (j + 1) * 4] = base_map.get(
                base, np.array([0, 0, 0, 0])
            )

    return matrix, descriptions


class Clusterer:
    def __init__(
        self,
        out_dir,
        reducer,
        clusterer,
        verbose=False,
        reducer_kwargs: dict = None,
        clusterer_kwargs: dict = None,
        encode: Literal["onehot", "pdist"] = "onehot",
    ):
        self.out_dir = out_dir
        self.reducer = reducer
        self.clusterer = clusterer
        self.add_config(verbose, reducer_kwargs, clusterer_kwargs, encode)

    def add_config(self, verbose, reducer_kwargs, clusterer_kwargs, encode):
        if encode not in ["onehot", "pdist"]:
            raise ValueError("encode must be one of 'onehot' or 'pdist'")
        self.config = config.Config()
        self.config.verbose = verbose

        if reducer_kwargs is None:
            reducer_kwargs = {}
        if clusterer_kwargs is None:
            clusterer_kwargs = {}
        self.config.add_seqcluster_config(
            reducer_kwargs, clusterer_kwargs, encode
        )

        fp_fh = base_logger.get_file_handler(
            os.path.join(self.out_dir, "cluster.log")
        )
        self.config.logger.addHandler(fp_fh)

    @base_logger.prog_log("Encode sequences")
    def encode_fasta(self, fa_path, encode):
        self.pdist_mx, self.descriptions = encode_pdist_mx(fa_path)
        if encode == "onehot":
            self.X, _ = encode_onehot_mx(fa_path)
        else:
            self.X = self.pdist_mx.copy()

    @base_logger.prog_log("Predict sequence clusters")
    def fit_predict(self, reducer_kwargs, clusterer_kwargs):
        if not hasattr(self.reducer, "fit_transform"):
            raise ValueError("Reducer does not have a fit_transform() method")
        if not hasattr(self.clusterer, "fit_predict"):
            raise ValueError("Clusterer does not have a fit_predict() method")
        self.embedding = self.reducer(**reducer_kwargs).fit_transform(self.X)
        self.cluster_labels = self.clusterer(**clusterer_kwargs).fit_predict(
            self.embedding
        )

    @base_logger.prog_log("Convert results to dataframe")
    def results_to_df(self):
        self.cluster_df = pd.DataFrame(
            {
                "description": self.descriptions,
                "cluster_label": self.cluster_labels,
                "embedding1": self.embedding[:, 0],
                "embedding2": self.embedding[:, 1],
            }
        )
        self.cluster_df.to_csv(
            os.path.join(self.out_dir, "cluster_df.csv"), index=False
        )

    def predict(self, data):
        if isinstance(data, BPData | MitoData):
            writer = write.Writer(data, self.config.n_cpu, self.config.verbose)
            with TemporaryDirectory() as tmp_dir:
                tmp_fasta = os.path.join(tmp_dir, "tmp.fasta")
                in_fasta = os.path.join(self.out_dir, "mltree.aln")
                writer.fasta(tmp_fasta)
                writer.align_fasta(tmp_fasta, in_fasta)
        elif isinstance(data, str):
            if not os.path.exists(data):
                raise FileNotFoundError(f"File not found: {data}")
            in_fasta = data
        else:
            raise ValueError(
                "data should be a BPData, MitoData class, or valid .FASTA file."
            )

        config = self.config.get_seqcluster_config()
        self.encode_fasta(in_fasta, config["encode"])
        self.fit_predict(config["reducer_kwargs"], config["clusterer_kwargs"])
        self.results_to_df()
        return self.cluster_df
