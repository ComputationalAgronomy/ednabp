from analysis_toolkit.runner_build import base_logger

class Reader:
    def __init__(self):
        pass


class BlastReader(Reader):
    DESIRED_LEVEL = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

    TAX_REPLACMENT = {"Mugil": "Mugilidae"}
                      #"KEY2": "REPLACE2"} # Allow multiple replacements later

    @staticmethod
    def generate_error_table(error_code:str = ':/\\*?"<>|', replace_symbol: str = '_') -> dict[str, str]:
        """
        Generate an error table to translate illegal characters in species names to a standard symbol.

        :param error_code: The error code used in the BLAST table.
        :param replace_symbol: The symbol to replace the error code with.
        :return: A dictionary containing the error code as the key and the replace symbol as the value.
        """
        key2 = list(error_code)
        error_translation = dict.fromkeys(key2, replace_symbol)
        error_symbol = str.maketrans(error_translation)
        return error_symbol

    def __init__(self):
        super().__init__()
        self.error_table = BlastReader.generate_error_table()
        self.hap2level = {}

    def process_line(self, line: str):
        """
        Process a single line from the BLAST CSV table.
        The line list should be in the format:
        0: Haplotype_id
        (not used)1: Subject accession
        2-8: species, genus, family, order, class, phylum, kingdom
        (not used)9-18: pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore (check specifiers: https://www.biostars.org/p/88944/#88949)

        :param line: The line to process.
        :return: A tuple containing the haplotype_id (e.g. "Zotu1") and a dictionary of taxonomic levels (e.g. {"species": "spcA", ...}).
        """
        line_list = line.split(',')

        level_list = [str(line_list[i]) for i in range(2, 9)]
        # identity = line_list[9]
        # length = line_list[14]
        # evalue = line_list[17]
        # bitscore = line_list[18]

        level_list[0] = level_list[0].translate(self.error_table)
        if level_list[2] in self.TAX_REPLACMENT:
            level_list[2] = self.TAX_REPLACMENT[level_list[2]]
        hap2level_entry = dict(zip(self.DESIRED_LEVEL, level_list))
        haplotype = line_list[0]
        self.update_hap2level(haplotype, hap2level_entry)

    def update_hap2level(self, haplotype, hap2level_entry):
        self.hap2level[haplotype] = hap2level_entry

    def read_blast_table(self, blast_table: str):
        """
        Read the BLAST CSV table and update the dictionary 'self.hap2level' with the corresponding taxonomic names at each level for every haplotype (ZOTU).
        Seven levels are used: species, genus, family, order, class, phylum, kingdom.

        :param blast_table_path: Path to the BLAST CSV table.
        """
        base_logger.logger.info(f"Reading Blast CSV Table: {blast_table}.")

        with open(blast_table, 'r') as file:
            for line in file.readlines():
                self.process_line(line)

        hap_count = len(self.hap2level)
        base_logger.logger.info(f"COMPLETE: Assigned {hap_count} haplotypes to species.")