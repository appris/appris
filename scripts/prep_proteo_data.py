#!/usr/bin/env python3
"""Prepare PROTEO data for use by APPRIS."""

from argparse import ArgumentParser
from collections import Counter
from csv import Dialect, DictReader, DictWriter, QUOTE_NONE, Sniffer
import os
import shutil
from tempfile import TemporaryDirectory


class ProteoCSV(Dialect):
    delimiter = ','
    doublequote = False
    escapechar = '\\'
    lineterminator = '\n'
    quotechar = '"'
    quoting = QUOTE_NONE
    skipinitialspace = False
    strict = True


ap = ArgumentParser(description=__doc__)
ap.add_argument('input_file',
                help='input proteomics data file')
ap.add_argument('output_file',
                help='prepared PROTEO CSV file')
args = ap.parse_args()

in_file_path = args.input_file
out_file_path = args.output_file


in_col_names = ('peptide', 'gene_id', 'mapped_transc_ids', 'unmapped_transc_ids',
                'best_pep_score', 'num_experiments')
out_col_names = ('peptide', 'gene_id', 'mapped_transc_ids', 'unmapped_transc_ids',
                 'num_experiments')

sniffer = Sniffer()
pep_seq_ctr = Counter()
known_delimiters = [',', '\t']
with TemporaryDirectory() as tmp_dir:
    tmp_file_path = os.path.join(tmp_dir, 'proteo.csv')

    with open(in_file_path) as in_f, open(tmp_file_path, 'w') as tmp_f:
        obs_dialect = sniffer.sniff(in_f.read(4096), known_delimiters)
        in_f.seek(0)
        reader = DictReader(in_f, in_col_names, dialect=obs_dialect)
        writer = DictWriter(tmp_f, out_col_names, extrasaction='ignore',
                            dialect=ProteoCSV)
        for row in reader:
            pep_seq_ctr.update([row['peptide']])
            writer.writerow(row)

    dup_pep_seqs = [x for x, n in pep_seq_ctr.items() if n > 1]

    if dup_pep_seqs:
        dup_pep_lines = '\n'.join(['  ' + x for x in dup_pep_seqs])
        raise ValueError(f"""duplicate peptide sequence(s):\n{dup_pep_lines}""")

    shutil.move(tmp_file_path, out_file_path)
