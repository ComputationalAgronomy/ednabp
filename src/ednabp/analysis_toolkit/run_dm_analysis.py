from .runner_exec.dm_barchart_plotter import BarchartPlotter
from .runner_exec.dm_contourf_plotter import ContourPlotter
from .runner_exec.dm_heatmap_plotter import HeatmapPlotter
from .runner_exec.dm_table_writer import DMWriter
from .runner_exec.dm_rankcorr_plotter import RankCorrPlotter

class DMAnalyser(BarchartPlotter, ContourPlotter, HeatmapPlotter, DMWriter, RankCorrPlotter):
    """
    A class that inherits from ContourPlotter, BarchartPlotter, HeatmapPlotter, and DMWriter.
    This class is used for diversity metrics related analysis and visualization.
    """
    def __init__(self, sampledata = None, no_verbose: bool = False):
        if sampledata is not None:
            DMWriter.__init__(self, sampledata, no_verbose)
        else:
            BarchartPlotter.__init__(self, no_verbose)