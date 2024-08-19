import os
import pandas as pd
import plotly.express as px

from analysis_toolkit.runner_build import (base_runner, utils)


class BarchartRunner(base_runner.AbundanceRunner):

    def __init__(self, samplesdata):
        super().__init__(samplesdata)

    def _fill_missing_keys(self):
        for sample_id in self.sample_id_used:
            self.samples2abundance[sample_id] = [self.samples2abundance[sample_id].get(unit_name, 0) for unit_name in self.uniq_unit_names]

    def _create_barchart_fig(self,
            axes_title_font: int = 20,
            axes_tick_font: int = 18,
            legend_font: int = 15,
            legend_x_position: float = 1.05,
            legend_y_position: float = 1.0
        ):
        """
        Create a stacked bar chart figure using Plotly.
        """
        plotdata = pd.DataFrame(self.samples2abundance, index=self.uniq_unit_names)
        data = plotdata.transpose()

        self.fig = px.bar(
            data,
            barmode='stack',
            labels={'value': 'Percentage (%)'},
            color_discrete_sequence=px.colors.qualitative.Pastel
        )
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

    def _save(self, save_html_dir: str, level):
        """
        Save the barchart as an HTML file.

        :param save_html_dir: The directory to save the HTML file.
        :param save_html_name: The name of the HTML file. If not provided, the name will be "{level}_barchart". Default is None.
        """
        bar_chart_path = os.path.join(save_html_dir, f"{level}_barchart.html")
        self.fig.write_html(bar_chart_path)
        self.logger.info(f"Barchart saved to: {bar_chart_path}")
 
        self.results_dir = save_html_dir

    def run_write(self):
        return super().run_write()

    @base_runner.log_execution("Plot relative abundance barchart", "plot_barchart.log")
    def run_plot(self,
            level: str,
            sample_id_list: list[str] = [],
            save_dir: str = None
        ):
        """
        Plot a barchart to visualize the abundance of a level across samples.

        :param level: The name of the level to plot (e.g., species, family, etc.).
        :param sample_id_list: A list of sample IDs to plot. Default is None (plot all samples).
        :param save_dir: If provided, the barchart will be saved as a .HTML file and save a log file. Default is None.
        """
        self._load_sample_id_list(sample_id_list)

        for sample_id in self.sample_id_used:
            self._load_units2abundance_dict(sample_id, level)
            self._normalize_abundance()
            self._update_samples2abundance_dict(sample_id)
        
        all_unit_names = [list(self.samples2abundance[sample_id].keys()) for sample_id in self.sample_id_used]
        self.uniq_unit_names = utils.list_union(all_unit_names)

        self._fill_missing_keys()
        
        self._create_barchart_fig()
        self.fig.show()

        if save_dir:
            self._save(save_dir, level)

        self.analysis_type = "barchart_plot"
        self.parameters.update(
            {
                "level": level,
            }
        )
    