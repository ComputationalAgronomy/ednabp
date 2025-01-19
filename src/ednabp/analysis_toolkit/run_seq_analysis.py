from .runner_exec.seq_mltree_writer import MLTreeWriter
from .runner_exec.seq_nexus_writer import NexusWriter
from .runner_exec.seq_umap_runner import UmapRunner

class SeqAnalyser(MLTreeWriter, NexusWriter, UmapRunner):
    pass