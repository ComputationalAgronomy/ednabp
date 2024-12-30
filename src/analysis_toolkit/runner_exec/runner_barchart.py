import os
import pandas as pd
import plotly.express as px

from analysis_toolkit.runner_build import base_runner


class BarchartRunner(base_runner.AbundanceRunner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)

    def run_write(self,
            write_type: str = "abundance",
            taxa_level: str = "species",
            save_dir: str = ".",
            normalize: bool = True,
            sample_id_list: list[str] = []
        ):
        return super().run_write(
            write_type=write_type,
            taxa_level=taxa_level,
            save_dir=save_dir,
            normalize=normalize,
            sample_id_list=sample_id_list
        )

    @base_runner.log_execution("Plot barchart", "plot_barchart.log")
    def run_plot(self,
            csv_path: str,
            save_dir: str = None,
            dereplicate: bool = False   
        ):
        """
        Plot a barchart to visualize the abundance of a level across samples.

        :param level: The name of the level to plot (e.g., species, family, etc.).
        :param sample_id_list: A list of sample IDs to plot. Default is None (plot all samples).
        :param save_dir: If provided, the barchart will be saved as a .HTML file and save a log file. Default is None.
        """
        self.df = pd.read_csv(csv_path)

        if dereplicate:
            self.df = self.df.drop("Sample", axis=1)

        self.df["Sample_id"] = self.df.apply(lambda x: "-".join(map(str, x[1:-1])), axis=1)

        self.fig = px.bar(self.df, x="Sample_id", y="Counts", color=self.df.columns[0])
        self._add_fig_setting()
        self.fig.show()

        if save_dir:
            self._save_html("Barchart", save_dir, os.path.basename(csv_path).split(".")[0])

        self.analysis_type = "Plot barchart"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "csv_path": csv_path,
                "dereplicate": dereplicate
            }
        )
    
    def _add_fig_setting(self,
            axes_title_font: int = 20,
            axes_tick_font: int = 18,
            legend_font: int = 15,
            legend_x_position: float = 1.05,
            legend_y_position: float = 1.0
        ):
        """
        Create a stacked bar chart figure using Plotly.
        """
        self.fig.update_xaxes(
            tickmode='linear',
            title=dict(
                text="Sample ID",
                font=dict(size=axes_title_font)
                ),
            tickfont=dict(size=axes_tick_font)
        )
        self.fig.update_yaxes(
            title=dict(
                text="Percentage (%)",
                font=dict(size=axes_title_font)
            ),
            tickfont=dict(size=axes_tick_font)
        )
        self.fig.update_layout(
            legend={
                "x": legend_x_position,
                "y": legend_y_position,
                "traceorder": 'normal',
                "orientation": 'h',
                "font": dict(size=legend_font)
            },
        )

        ## https://stackoverflow.com/questions/44309507/stacked-bar-plot-using-matplotlib
        # a way to sort stacked BarChart
        # plotdata = plotdata.transpose()
        # fig, ax = plt.subplots()
        # x = plotdata.index
        # indexes = np.argsort(plotdata.values).T
        # heights = np.sort(plotdata.values).T
        # order = -1
        # bottoms = heights[::order].cumsum(axis=0)
        # bottoms = np.insert(bottoms, 0, np.zeros(len(bottoms[0])), axis=0)
        # colormap = plt.get_cmap('rainbow')
        # colors = px.colors.qualitative.Pastel
        # num_colors = len(colors)
        # mpp_colors = {col: colors[i % num_colors] for i, col in enumerate(plotdata.columns)}
    
        # fig = go.Figure()
        
        # # # Plot each segment of the stacked bar
        # for i, (idxs, vals) in enumerate(list(zip(indexes, heights))[::order]):
        #     mps = np.take(np.array(plotdata.columns), idxs)
        #     for j, m in enumerate(mps):
    
        #         fig.add_trace(go.Bar(
        #             x=x,
        #             y=[vals[j]] * len(x),
        #             name=m,
        #             marker=dict(color=mpp_colors[m]),
        #             offsetgroup=i,
        #             base=bottoms[i][j]
        #         ))
    
        # # Update layout
        # fig.update_layout(
        #     barmode='stack',
        #     xaxis_title="Sample ID",
        #     yaxis_title="Percentage (%)",
        #     legend=dict(
        #         x=1.05,
        #         y=1,
        #         traceorder='normal',
        #         orientation='h',
        #         font=dict(size=15)
        #     ),
        # )
        # fig.update_xaxes(tickmode='linear', title=dict(font=dict(size=20)), tickfont=dict(size=18))
        # fig.update_yaxes(title=dict(font=dict(size=20)), tickfont=dict(size=18))
        # fig.show()