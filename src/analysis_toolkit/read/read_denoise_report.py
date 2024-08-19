import re

from analysis_toolkit.read import read_blast_csv
from analysis_toolkit.runner_build import base_logger

class DenoiseReportReader(read_blast_csv.Reader):
    RE_DENOISE_PATTERN = re.compile(r"(Uniq\d*)|size=(\d*)|(amp\d*)|top=(Uniq\d*)")
    RE_CHFILTER_PATTERN = re.compile(r"(Uniq\d*)|size=(\d*)|(zotu)|(chimera)")

    def __init__(self):
        super().__init__()
        self.amp_size = {}
        self.hap2amp = {}
        self.hap_size = {}

    @staticmethod
    def read_denoise_line(line: str) -> list[str]:
        """
        Read a line containing 'denoise' string and return a list of the values.

        :param line: one line from the denoise report
        :return: [Amplicon_id, Size, Top]. e.g. ['Uniq1', '88422', 'amp1'], ['Uniq6', '8126', 'Uniq2']
        """
        line_list = [
            "".join(t)
            for t in DenoiseReportReader.RE_DENOISE_PATTERN.findall(line)
        ]

        # If Top is "ampXX", it is a true amplicon, otherwise, it is a noise.
        if 'amp' in line_list[2]:
            line_list[2] = line_list[0]

        return line_list

    @staticmethod
    def read_chifilter_line(line: str) -> list[str]:
        """
        Read a line containing 'chifilter' string and return a list of the values.

        :param line: one line from the denoise report
        :return: [Amplicon_id, Size, Assigned_type]. e.g. ['Uniq1', '88422', 'zotu'], ['Uniq102', '124', 'chimera']
        """
        line_list = [
            "".join(t)
            for t in DenoiseReportReader.RE_CHFILTER_PATTERN.findall(line)
        ]

        return line_list
    
    def process_denoise_line(self, line: str):
        """
        Process a line containing 'denoise' string and update amp_size and hap2amp dictionaries.
        This step create relationship between haplotypes(top) and amplicons and record the size of each amplicon (it means unique sequence here).
        
        :param line: One line from the denoise report
        """
        amplicon, size, top = DenoiseReportReader.read_denoise_line(line)

        self.amp_size[amplicon] = size

        if top not in self.hap2amp:
            self.hap2amp[top] = []
        self.hap2amp[top].append(amplicon)

    def process_chifilter_line(self, line: str, zotu_count: int, chimera_count: int) -> tuple[int, int]:
        """
        Process a line containing 'chifilter' string and update hap_size and hap2amp dictionaries.
        This step rename the top from 'UniqXX' to 'ZotuYY' or 'ChimeraZZ' and record the size for each haplotype (it means ZOTU here).
        
        :param line: One line from the denoise report
        :param zotu_count: Current count of zotus
        :param chimera_count: Current count of chimeras
        :return: Updated counts of zotus and chimeras
        """
        old_top, size, assigned_type = DenoiseReportReader.read_chifilter_line(line)

        if assigned_type == 'zotu':
            zotu_count += 1
            new_top = f'Zotu{zotu_count}'
            self.hap_size[new_top] = size
        elif assigned_type == 'chimera':
            chimera_count += 1
            new_top = f'Chimera{chimera_count}'
        
        if old_top in self.hap2amp:
            self.hap2amp[new_top] = self.hap2amp.pop(old_top)
        
        return zotu_count, chimera_count

    def read_denoise_report(self, denoise_report: str):
        """
        read a denoise report (.txt) generated from usearch and update amp_size, hap2amp and hap_size dictionaries.

        :param denoise_report_path: path to the denoise report
        """
        base_logger.logger.info(f"Reading denoising Report:  {denoise_report}")

        zotu_count = 0
        chimera_count = 0

        with open(denoise_report, 'r') as file:
            for line in file.readlines():
                if 'denoise' in line:
                    self.process_denoise_line(line)
                elif 'chfilter' in line:
                    zotu_count, chimera_count = self.process_chifilter_line(line, zotu_count, chimera_count)

        base_logger.logger.info(f"COMPLETE: Total {zotu_count} biological Haplotypes kept; Total {chimera_count} predicted Chimeras removed.")