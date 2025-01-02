import os
import pandas as pd
import plotly.express as px

from ..runner_build import DMRunner, base_logger


class BarchartRunner(DMRunner):

    def __init__(self, sampledata, no_verbose):
        super().__init__(sampledata, no_verbose)

    @base_logger.prog_log("Plot barchart")
    def plot_barchart(self,
            csv_path: str,
            taxa_column: str,
            metric_column: str,
            save_dir: str | None = None,
            overwrite: bool = False
        ):
        """
        Plot a barchart to visualize the abundance of a level across samples.

        :param csv_path: Path to the CSV file containing the data
        :param taxa_column: Column name to use for color values 
        :param metric_column: Column name to use for y-axis values
        :param save_dir: If provided, the barchart will be saved as a .HTML file. Default is None.
        :param overwrite: If True, overwrites existing files in save_dir. If False, don't overwrite. Defaults to False.
        """
        self._load_and_validate_data(csv_path, taxa_column, metric_column)
        self._create_pivot_table(taxa_column, metric_column)
        self._prepare_plot_data()
        self._create_plot()
        self._display_and_save_plot(save_dir, overwrite)

    @base_logger.prog_log("Load and validate data")
    def _load_and_validate_data(self, csv_path: str, taxa_column: str, metric_column: str):
        try:
            self.df = pd.read_csv(csv_path)
            required_columns = {'Sample_id', taxa_column, metric_column}
            if not required_columns.issubset(self.df.columns):
                missing = required_columns - set(self.df.columns)
                raise ValueError(f"Missing required columns: {missing}")
        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find CSV file: {csv_path}")

    @base_logger.prog_log("Process data")
    def _create_pivot_table(self, taxa_column: str, metric_column: str):
        self.pivot_df = self.df.pivot(
            index=self.SAMPLE_ID_COLUMN,
            columns=taxa_column,
            values=metric_column
        ).fillna(0)

        # Sort columns by sum of values (descending)
        column_sums = self.pivot_df.sum()
        self.pivot_df = self.pivot_df[column_sums.sort_values(ascending=False).index]

        # Sort rows by sum of values (descending)
        row_sums = self.pivot_df.sum(axis=1)
        self.pivot_df = self.pivot_df.loc[row_sums.sort_values(ascending=False).index]

    @base_logger.prog_log("Prepare plot data")
    def _prepare_plot_data(self):
        cols = self.pivot_df.columns.tolist()
        idx = self.pivot_df.index.tolist()
        self.plot_data = {
            'x': [i for i in idx for _ in range(len(cols))],
            'y': self.pivot_df.values.flatten(),
            'color': cols * len(idx)
        }

    @base_logger.prog_log("Create plot")
    def _create_plot(self):
        self.fig = px.bar(
            x=self.plot_data['x'],
            y=self.plot_data['y'],
            color=self.plot_data['color']
        )
        self._add_fig_setting()

    def _add_fig_setting(self,
            axes_title_font: int = 20,
            axes_tick_font: int = 18,
            legend_font: int = 15,
            legend_x_position: float = 1.05,
            legend_y_position: float = 1.0
        ):
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

    @base_logger.prog_log("Display and save plot (if 'save_dir' provided")
    def _display_and_save_plot(self, save_dir: str | None, overwrite: bool):
        self.fig.show()
        if save_dir:
            self._save_html(save_dir, overwrite)

    def _save_html(self, save_html_dir: str, overwrite: bool):
        fig_path = os.path.join(save_html_dir, "barchart.html")
        if os.path.exists(fig_path) and not overwrite:
            self.logger.warning(f"WARNING: File already exists: {fig_path}. Stop saving.")
            return
        self.fig.write_html(fig_path)
        self.logger.info(f"Barchart saved to: {fig_path}")