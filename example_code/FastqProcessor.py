# cli: (eDNAbp -h)
# example: eDNAbp -i .\stage_test\fastq\ -o .\stage_test\ -db .\data\database\MiFish -ln .\data\database\lineage.csv --n-cpu 20
from ednabp import FastqProcessor

FastqProcessor(
    input_path=".\\stage_test\\fastq",
    output_path=".\\stage_test",
    db_path=".\\data\\database\\MiFish",
    lineage_path=".\\data\\database\\lineage.csv",
    n_cpu=20,
)