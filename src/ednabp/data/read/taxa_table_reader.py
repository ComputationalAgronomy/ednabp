import pandas as pd

from ...common import base_logger
from .base_reader import Reader


def generate_error_table(
    error_code: str = ':/\\*?"<>|', replace_symbol: str = "_"
) -> dict[str, str]:
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


class TaxaTableReader(Reader):
    DESIRED_LEVEL = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
    ]

    TAX_REPLACMENT = {"Mugil": "Mugilidae"}
    # "KEY2": "REPLACE2"} # Allow multiple replacements later

    def __init__(self):
        super().__init__()
        self.error_table = generate_error_table()
        self.hap2level = {}

    @base_logger.prog_log(prog_name="Read Taxa CSV Table")
    def read_taxa_table(self, taxa_table: str):
        df = pd.read_csv(taxa_table)
        key_column = "qseqid"
        value_columns = self.DESIRED_LEVEL
        df_indexed = df.set_index(key_column)[value_columns]
        self.hap2level = df_indexed.to_dict("index")
