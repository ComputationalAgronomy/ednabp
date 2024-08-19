import os
import subprocess
import tempfile

from analysis_toolkit.runner_build import (base_logger, utils)


def derep_fasta(seq_path: str, uniq_path: str, relabel: str, threads: int = 12, sizeout: bool = False) -> None:
    """
    Dereplicate a FASTA file by removing duplicate sequences.
    (USEARCH command reference: https://drive5.com/usearch/manual/cmd_fastx_uniques.html)

    :param seq_path: Path to the input FASTA file.
    :param uniq_path: Path to the output FASTA file with unique sequences.
    :param relabel: Prefix to add to the sequence labels.
    :param threads: Number of threads to use. Default is 12.
    :param sizeout: If True, size annotations will be added to the output sequence labels. Default is False.
    """
    cmd = [
        'usearch', '-fastx_uniques', seq_path, '-threads', str(threads),
        '-relabel', f'{relabel}-', '-fastaout', uniq_path
    ]
    if sizeout:
        cmd.append('-sizeout')

    utils.run_subprocess("USEARCH", cmd, uniq_path)

def write_fasta(units2fasta_dict: dict[str, str], save_path: str, dereplicate: bool = False, sizeout: bool = False) -> None:
    """
    Write sequences to a FASTA file.

    :param units2fasta_dict: Dictionary with sequence names (unit names) as keys and FASTA format sequences as values.
    :param save_path: Path to the output FASTA file.
    :param dereplicate: If True, dereplicate the sequences before writing to the file. Default is False.
    :param sizeout: If True, size annotations will be added to the output sequence labels. Only works when dereplicating. Default is False.
    """
    if dereplicate:
        fasta_list = []
        with tempfile.TemporaryDirectory() as tmpdirname:
            for unit_name, unit_seq in units2fasta_dict.items():
                unit_fa_path = os.path.join(tmpdirname, f'{unit_name}.fa')
                unit_uniq_fa_path = os.path.join(tmpdirname, f'{unit_name}_uniq.fa')

                with open(unit_fa_path, 'w') as file:
                    file.write(unit_seq)

                derep_fasta(seq_path=unit_fa_path, uniq_path=unit_uniq_fa_path, relabel=unit_name, sizeout=sizeout)

                with open(unit_uniq_fa_path, 'r') as file:
                    fasta_list.append(file.read())
    else:
        fasta_list = list(units2fasta_dict.values())

    fasta_str = "".join(fasta_list)
    with open(save_path, 'w') as file:
        file.write(fasta_str)

    num_seq = fasta_str.count(">")
    base_logger.logger.info(f"Written {num_seq} sequences to: {save_path}")

    return num_seq

def align_fasta(seq_path: str, aln_path: str) -> None:
    """
    Align sequences in a FASTA file using Clustal Omega and save the aligned sequences to an output file.
    (ClustalO command reference: http://www.clustal.org/omega/README)

    :param seq_path: Path to the input FASTA file containing sequences to align.
    :param aln_path: Path to the output FASTA file to save the aligned sequences.
    """
    cmd = [
        'clustalo', '-i', seq_path, '-o', aln_path, '--force'
    ]

    utils.run_subprocess("Clustal Omega", cmd, aln_path)