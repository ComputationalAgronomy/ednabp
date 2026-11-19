import os
import shutil
import subprocess

from ...common.default_settings import SETTINGS
from ..step_build import stage_builder
from .otutab import add_otutab_stages
from .utils import read_fasta, write_fasta


def parse_uc(uc_file):
    clusters = {}
    with open(uc_file) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                clusters.setdefault(int(fields[1]), []).insert(0, fields[8])
            elif fields[0] == "H":
                clusters.setdefault(int(fields[1]), []).append(fields[8])
    return clusters


def size_of(label):
    fields = dict(p.split("=", 1) for p in label.split(";") if "=" in p)
    return int(fields["size"])


def parse_consensus_labels(consout_file):
    drafts = {}
    for label, seq in read_fasta(consout_file):
        centroid = label.removeprefix("centroid=")
        end = centroid.rfind(";seqs=")
        drafts[centroid[:end] if end >= 0 else centroid] = seq
    return drafts


class ClusterOntStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        vsearch_prog="vsearch",
        minimap2_prog="minimap2",
        racon_prog="racon",
        medaka_prog="medaka_consensus",
        usearch_prog="usearch",
        in_dir="",
        out_dir="",
        in_suffix=".fasta",
        out_suffix=".fasta",
        cluster_id: float = 0.90,
        min_cluster_size: int = 10,
        max_reads: int = 200,
        racon_rounds: int = 2,
        medaka: bool = False,
        medaka_model: str | None = None,
        consensus_id: float = 0.99,
        cluster_strand: bool = False,
        otutab_id: float = 0.90,
        otutab_strand: bool = True,
        keep_tmp: bool = False,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.vsearch_prog = vsearch_prog
        self.minimap2_prog = minimap2_prog
        self.racon_prog = racon_prog
        self.medaka_prog = medaka_prog
        self.usearch_prog = usearch_prog
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.cluster_id = cluster_id
        self.min_cluster_size = min_cluster_size
        self.max_reads = max_reads
        self.racon_rounds = racon_rounds
        self.medaka = medaka
        self.medaka_model = medaka_model
        self.consensus_id = consensus_id
        self.cluster_strand = cluster_strand
        self.otutab_id = otutab_id
        self.otutab_strand = otutab_strand
        self.keep_tmp = keep_tmp
        self.summary_suffix = SETTINGS["suffix"]["cluster_ont_summary"]
        self.report_suffix = SETTINGS["suffix"]["report"]

    def setup(self, prefix):
        self.prefix = prefix
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.summary_file = os.path.join(
            self.out_dir, f"{prefix}{self.summary_suffix}"
        )
        self.tmp_dir = os.path.join(self.out_dir, f"{prefix}_cluster_ont_tmp")
        self.uc_file = os.path.join(self.tmp_dir, "clusters.uc")
        self.consout_file = os.path.join(self.tmp_dir, "draft_consensus.fasta")
        report = os.path.join(self.out_dir, f"{prefix}{self.report_suffix}")
        self.check_infile()
        os.makedirs(self.tmp_dir, exist_ok=True)

        strand = "both" if self.cluster_strand else "plus"
        cmd = (
            f"{self.vsearch_prog} --cluster_fast {self.infile}"
            f" --id {self.cluster_id} --strand {strand}"
            f" --threads {self.config.n_cpu}"
            f" --uc {self.uc_file} --consout {self.consout_file}"
        )
        super().add_stage("Cluster ONT reads", cmd)
        super().add_stage_output_to_file(
            "Write vsearch report", 0, report, report
        )
        super().add_stage_function(
            "Polish cluster consensus sequences", self.polish_clusters
        )
        add_otutab_stages(
            self,
            prefix,
            self.infile,
            self.outfile,
            zotus=False,
            identity=self.otutab_id,
            strand=self.otutab_strand,
        )
        if not self.keep_tmp:
            super().add_stage_function(
                "Remove ONT clustering temporary files",
                lambda: shutil.rmtree(self.tmp_dir, ignore_errors=True),
            )

    def run_cmd(self, args, stdout_path=None):
        if self.config.verbose:
            self.config.logger.info(f"Command: {' '.join(args)}")
        if stdout_path is None:
            subprocess.run(args, check=True, capture_output=True, text=True)
        else:
            with open(stdout_path, "w") as out:
                subprocess.run(
                    args,
                    check=True,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                )

    def racon_polish(self, reads_file, draft_file, work_dir):
        for i in range(1, self.racon_rounds + 1):
            paf = os.path.join(work_dir, f"round{i}.paf")
            polished = os.path.join(work_dir, f"racon{i}.fasta")
            self.run_cmd(
                [
                    self.minimap2_prog,
                    "-x",
                    "map-ont",
                    "-t",
                    str(self.config.n_cpu),
                    draft_file,
                    reads_file,
                ],
                stdout_path=paf,
            )
            if os.path.getsize(paf) == 0:
                break  # no read aligned; keep the current draft
            self.run_cmd(
                [
                    self.racon_prog,
                    "-t",
                    str(self.config.n_cpu),
                    reads_file,
                    paf,
                    draft_file,
                ],
                stdout_path=polished,
            )
            if os.path.getsize(polished) == 0:
                break
            draft_file = polished
        return draft_file

    def medaka_polish(self, reads_file, draft_file, work_dir):
        medaka_dir = os.path.join(work_dir, "medaka")
        args = [
            self.medaka_prog,
            "-i",
            reads_file,
            "-d",
            draft_file,
            "-o",
            medaka_dir,
            "-t",
            str(self.config.n_cpu),
        ]
        if self.medaka_model:
            args += ["-m", self.medaka_model]
        self.run_cmd(args)
        return os.path.join(medaka_dir, "consensus.fasta")

    def polish_clusters(self):
        reads = {
            label.split()[0]: seq for label, seq in read_fasta(self.infile)
        }
        clusters = parse_uc(self.uc_file)
        drafts = parse_consensus_labels(self.consout_file)

        kept = sorted(
            (
                members
                for members in clusters.values()
                if len(members) >= self.min_cluster_size
            ),
            key=len,
            reverse=True,
        )
        n_reads = sum(len(m) for m in clusters.values())
        self.config.logger.info(
            f"{len(clusters)} clusters from {n_reads} reads;"
            f" {len(kept)} clusters with >= {self.min_cluster_size} reads"
            f" ({sum(len(m) for m in kept)} reads) will be polished."
        )

        consensus = []
        for i, members in enumerate(kept, 1):
            work_dir = os.path.join(self.tmp_dir, f"cluster{i}")
            os.makedirs(work_dir, exist_ok=True)
            reads_file = os.path.join(work_dir, "reads.fasta")
            draft_file = os.path.join(work_dir, "draft.fasta")
            write_fasta(
                ((r, reads[r]) for r in members[: self.max_reads]), reads_file
            )
            draft_seq = drafts.get(members[0], reads[members[0]])
            write_fasta([(f"cluster{i}", draft_seq)], draft_file)

            polished = self.racon_polish(reads_file, draft_file, work_dir)
            if self.medaka:
                polished = self.medaka_polish(reads_file, polished, work_dir)

            seq = next(read_fasta(polished))[1]
            consensus.append((f"cluster{i};size={len(members)}", seq))

        if not consensus:
            raise RuntimeError(
                f"No cluster has >= {self.min_cluster_size} reads;"
                " lower min_cluster_size or check the input reads."
            )
        self.merge_consensus(consensus)

    def merge_consensus(self, consensus):
        """Merge near-identical consensuses and write the final OTUs."""
        consensus_file = os.path.join(self.tmp_dir, "polished.fasta")
        centroids_file = os.path.join(self.tmp_dir, "merged.fasta")
        merge_uc = os.path.join(self.tmp_dir, "merged.uc")
        write_fasta(consensus, consensus_file)
        self.run_cmd(
            [
                self.vsearch_prog,
                "--cluster_size",
                consensus_file,
                "--id",
                str(self.consensus_id),
                "--strand",
                "both" if self.cluster_strand else "plus",
                "--sizein",
                "--sizeout",
                "--centroids",
                centroids_file,
                "--uc",
                merge_uc,
                "--threads",
                str(self.config.n_cpu),
            ]
        )

        merged = {
            labels[0].split(";")[0]: [m.split(";")[0] for m in labels]
            for labels in parse_uc(merge_uc).values()
        }
        otus = sorted(
            (
                (size_of(label), label.split(";")[0], seq)
                for label, seq in read_fasta(centroids_file)
            ),
            key=lambda x: x[0],
            reverse=True,
        )

        with (
            open(self.outfile, "w") as fa,
            open(self.summary_file, "w") as tsv,
        ):
            tsv.write("otu\tcluster_reads\tlength\tmerged_clusters\n")
            for i, (size, cluster, seq) in enumerate(otus, 1):
                fa.write(f">Otu{i}\n{seq}\n")
                tsv.write(
                    f"Otu{i}\t{size}\t{len(seq)}\t{','.join(merged[cluster])}\n"
                )
        self.config.logger.info(
            f"{len(consensus)} polished consensus sequences merged into"
            f" {len(otus)} OTUs at {self.consensus_id} identity."
        )

    def run(self):
        super().run()
        return all(self.output)


def cluster_ont_demo(config, prefix, reads_dir="", save_dir=""):
    stage = ClusterOntStage(config, in_dir=reads_dir, out_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
