#!/usr/bin/env python3
"""Restore Gencode/Ensembl ID versions in TSV file."""

from argparse import ArgumentParser
from contextlib import contextmanager
import csv
import functools
import gzip
import os
import re
from string import ascii_uppercase, digits

from gtfparse import read_gtf


ensembl_id_regex = re.compile(
    '^(?P<id>ENS(?:[A-Z]{3})?(?:G|T|E|P)\d{11})'
    '(?P<version>\.\d+(?P<par_y_suffix>_PAR_Y)?)?$'
)


@contextmanager
def open_as_text(file, mode='r'):
    if mode not in {'a', 'r', 'w', 'x'}:
        raise ValueError(f"unknown file mode: '{mode}'")
    text_mode = mode + 't'
    if str(file).lower().endswith('.gz'):
        opener = functools.partial(gzip.open, mode=text_mode)
    else:
        opener = functools.partial(open, mode=text_mode)
    file_obj = None
    try:
        file_obj = opener(file)
        yield file_obj
    finally:
        if file_obj is not None:
            file_obj.close()


def strip_ensembl_id_version(ensembl_id):
    match = ensembl_id_regex.match(ensembl_id)
    try:
        bare_id = match['id']
    except TypeError:
        raise ValueError(f"unrecognised Ensembl ID: '{ensembl_id}'")
    return bare_id


if __name__ == '__main__':

    ap = ArgumentParser(description=__doc__)
    ap.add_argument('-i', '--input-file', metavar='PATH', required=True,
                    help='input TSV file')
    ap.add_argument('--gtf-file', metavar='PATH', required=True,
                    help='Gencode annotation file')
    ap.add_argument('--subfield-sep', metavar='STRING',
                    help='subfield separator, if a field can contain subfields')
    ap.add_argument('-o', '--output-file', metavar='PATH', required=True,
                    help='output TSV file with ID versions')
    args = ap.parse_args()

    input_file = args.input_file
    gtf_file = args.gtf_file
    subfield_sep = args.subfield_sep
    output_file = args.output_file

    if subfield_sep is not None and (len(subfield_sep) != 1 or
            subfield_sep in ascii_uppercase + digits):
        raise ValueError(f"invalid subfield separator: '{subfield_sep}'")

    df = read_gtf(gtf_file)
    id_ver_map = dict()
    for k in ('gene_id', 'transcript_id', 'protein_id', 'exon_id'):
        id_ver_map.update({
            strip_ensembl_id_version(x): x for x in df[k]
            if x != '' and not x.endswith('_PAR_Y')
        })

    with open_as_text(input_file) as in_f, open_as_text(output_file, mode='w') as out_f:
        reader = csv.reader(in_f, dialect='excel-tab')
        writer = csv.writer(out_f, dialect='excel-tab',
                            lineterminator=os.linesep)
        for row in reader:
            for i, field in enumerate(row):
                match = ensembl_id_regex.match(field)
                if match:

                    id_unknown = False
                    if match['version'] is None:
                        try:
                            row[i] = id_ver_map[field]
                        except KeyError:
                            id_unknown = True
                    elif match['id'] not in id_ver_map:
                        id_unknown = True

                    if id_unknown:
                        raise ValueError(
                            f"ID '{field}' not found in provided Gencode annotation")

                elif subfield_sep is not None:
                    subfields = field.split(subfield_sep)
                    num_changes = 0
                    for j, subfield in enumerate(subfields):
                        match = ensembl_id_regex.match(subfield)

                        id_unknown = False
                        if match:
                            if match['version'] is None:
                                try:
                                    subfields[j] = id_ver_map[subfield]
                                except KeyError:
                                    id_unknown = True
                                else:
                                    num_changes += 1
                            elif match['id'] not in id_ver_map:
                                id_unknown = True

                            if id_unknown:
                                raise ValueError(
                                    f"ID '{subfield}' not found in provided Gencode annotation")

                    if num_changes > 0:
                        row[i] = subfield_sep.join(subfields)
            writer.writerow(row)
