import re

from .base_reader import Reader
from ..runner_build import base_logger


class DenoiseReportReader(Reader):
    RE_DENOISE_PATTERN = re.compile(r"(Uniq\d*)|size=(\d*)|(amp\d*)|top=(Uniq\d*)")
    RE_CHFILTER_PATTERN = re.compile(r"(Uniq\d*)|size=(\d*)|(zotu)|(chimera)")

    def __init__(self):
        super().__init__()
        self.amp_size = {}
        self.hap2amp = {}
        self.hap_size = {}
        self.zotu_count = 1
        self.chi_count = 1

    @base_logger.prog_log(prog_name="Read Denoise Report")
    def read_denoise_report(self, denoise_report: str):
        """
        read a denoise report (.txt) generated from usearch and update amp_size, hap2amp and hap_size dictionaries.

        :param denoise_report_path: path to the denoise report
        """
        with open(denoise_report, 'r') as file:
            for line in file.readlines():
                if 'denoise' in line:
                    self.process_denoise_line(line)
                elif 'chfilter' in line:
                    self.process_chifilter_line(line)

    def read_denoise_line(line: str) -> list[str]:
        """
        Read a line containing 'denoise' string and return a list of the values.

        :param line: one line from the denoise report
        :return: [Amplicon_id, Size, Top]. e.g. ['Uniq1', '88422', 'amp1'], ['Uniq6', '8126', 'Uniq2']
        """
        line_list = ["".join(t) for t in DenoiseReportReader.RE_DENOISE_PATTERN.findall(line)]
        # If Top is "ampXX", it is a true amplicon, otherwise, it is a noise.
        if 'amp' in line_list[2]:
            line_list[2] = line_list[0]
        return line_list

    def read_chifilter_line(line: str) -> list[str]:
        """
        Read a line containing 'chifilter' string and return a list of the values.

        :param line: one line from the denoise report
        :return: [Amplicon_id, Size, Assigned_type]. e.g. ['Uniq1', '88422', 'zotu'], ['Uniq102', '124', 'chimera']
        """
        line_list = ["".join(t) for t in DenoiseReportReader.RE_CHFILTER_PATTERN.findall(line)]
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

    def process_chifilter_line(self, line: str):
        """
        Process a line containing 'chifilter' string and update hap_size and hap2amp dictionaries.
        This step rename the top from 'UniqXX' to 'ZotuYY' or 'ChimeraZZ' and record the size for each haplotype (it means ZOTU here).
        
        :param line: One line from the denoise report
        """
        old_top, size, assigned_type = DenoiseReportReader.read_chifilter_line(line)

        if assigned_type == 'zotu':
            new_top = f'Zotu{self.zotu_count}'
            self.hap_size[new_top] = size
            self.zotu_count += 1
        elif assigned_type == 'chimera':
            new_top = f'Chimera{self.chi_count}'
            self.chi_count += 1
        else:
            raise ValueError(f"Unknown assigned_type: {assigned_type}")
        
        if old_top in self.hap2amp:
            self.hap2amp[new_top] = self.hap2amp.pop(old_top)