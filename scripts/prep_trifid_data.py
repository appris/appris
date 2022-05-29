#!/usr/bin/env python3
"""Prepare TRIFID prediction file for APPRIS.

This script fetches TRIFID prediction data from
either a local file, or if that is not specified,
the TRIFID GitLab repository, then sorts,
compresses and indexes the data.

Rows are sorted first
by gene ID, then in order of decreasing TRIFID score,
and then by transcript ID to break ties. Finally, a
transcript index column is added.

The file is then compressed with
bgzip and indexed by tabix.
"""

from argparse import ArgumentParser
from collections import defaultdict, namedtuple, OrderedDict
from collections.abc import MutableMapping
from contextlib import contextmanager
from datetime import datetime
from decimal import Decimal
import csv
import functools
import gzip
import hashlib
import os
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory

from Bio.SeqIO.FastaIO import SimpleFastaParser
from dateutil import parser
#from gitlab import Gitlab
from gtfparse import read_gtf
from pytz import timezone


class KeyClashError(KeyError):
    pass


class SlottedDict(MutableMapping):

    def __init__(self, *args, **kwargs):
        self._dict = dict()
        self.update(*args, **kwargs)

    def __delitem__(self, key):
        del self._dict[key]

    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        for key in self._dict:
            yield key

    def __len__(self):
        return len(self._dict)

    def __repr__(self):
        contents = repr(self._dict) if len(self._dict) > 0 else ''
        return f'{type(self).__name__}({contents})'

    def __setitem__(self, key, value):
        if key in self._dict:
            raise KeyClashError(repr(key))
        self._dict[key] = value

    def __str__(self):
        contents = self._dict if len(self._dict) > 0 else ''
        return f'{type(self).__name__}({contents})'

    def items(self):
        for item in self._dict.items():
            yield item

    def keys(self):
        for key in self._dict.keys():
            yield key

    def to_dict(self):
        return self._dict.copy()

    def values(self):
        for value in self._dict.values():
            yield value


def calc_sha1(s):
    return hashlib.sha1(s.encode('utf-8')).hexdigest()


def copy_local_trifid_file(local_file_path, pred_file_path):
    mod_dt = datetime.fromtimestamp(os.path.getmtime(local_file_path))
    pred_file_version = mod_dt.strftime('%Y%m%d')
    shutil.copyfile(local_file_path, pred_file_path)
    return pred_file_version


def download_gitlab_trifid_file(assembly, source_slug, pred_file_path):

    gl = Gitlab('https://gitlab.com')
    project = gl.projects.get('bu_cnio/trifid')

    # TODO: use pagination, allow date ranges
    repo_file_path = f'data/genomes/{assembly}/{source_slug}/trifid_predictions.tsv.gz'
    print(f"Fetching commit history of TRIFID file: '{repo_file_path}' …")
    query_params = {'ref_name': 'master', 'path': repo_file_path}
    commits = project.commits.list(all=True, query_parameters=query_params)

    if not commits:
        raise ValueError('no TRIFID prediction file found with given parameters')

    print('Selecting most recent version of TRIFID file …')
    latest_commit = None
    latest_commit_date = None
    for commit in commits:
        commit_date = get_utc_commit_date(commit)
        if latest_commit is None or commit_date > latest_commit_date:
            # TODO: implement option to download earlier version of prediction file
            latest_commit_date = commit_date
            latest_commit = commit
        elif commit_date == latest_commit_date:
            # TODO: handle same-day commits
            raise ValueError('cannot resolve unique latest commit date of TRIFID prediction file')

    pred_file_meta = project.files.get(repo_file_path, latest_commit.id)
    pred_file_version = latest_commit_date.strftime('%Y%m%d')

    print('Downloading TRIFID file …')
    with open(pred_file_path, 'wb') as f:
        project.repository_raw_blob(pred_file_meta.blob_id, streamed=True, action=f.write)

    return pred_file_version


def get_utc_commit_date(commit):
    ts_local = commit.committed_date
    dt_local = parser.parse(ts_local)
    dt_utc = dt_local.astimezone(timezone('UTC'))
    return dt_utc.date()


def is_par_y(gencode_id):
    return gencode_id.endswith('_PAR_Y')


def load_ensembl_metadata(annot_file, transl_file):

    transl_seq_meta = SlottedDict()
    with open_as_text(transl_file) as f:
        for header, sequence in SimpleFastaParser(f):
            match = ensembl_transl_hdr_regex.match(header)
            transl_id_ver = match['ens_transl_id']

            match = ens_id_regex.match(transl_id_ver)
            transl_key = (match['bare_id'], match['id_version'])

            transl_seq_meta[transl_key] = {
                'translation_seq_sha1': calc_sha1(sequence),
                'length': len(sequence)
            }

    transl_keys = set(transl_seq_meta.keys())

    annot_df = read_gtf(annot_file)

    annot_df = annot_df.loc[(annot_df['feature'] == 'CDS') &
                            (annot_df['protein_id'] != '') &
                            (annot_df['protein_version'] != '')]

    rel_cols = ['gene_id', 'gene_version', 'transcript_id',
                'transcript_version', 'protein_id', 'protein_version']
    annot_df = annot_df.loc[:, rel_cols].drop_duplicates(ignore_index=True)

    gene_to_isoforms = defaultdict(set)
    for _,row in annot_df.iterrows():
        gene_key = (row['gene_id'], row['gene_version'])
        transc_key = (row['transcript_id'], row['transcript_version'])
        transl_key = (row['protein_id'], row['protein_version'])
        gene_to_isoforms[gene_key].add((transc_key, transl_key))

    isoform_meta = SlottedDict()
    for gene_key, isoforms in gene_to_isoforms.items():
        gene_id, gene_version = gene_key

        for transc_key, transl_key in isoforms:
            transc_id, transc_version = transc_key
            transl_id, transl_version = transl_key

            transl_seq_sha1 = transl_seq_meta[transl_key]['translation_seq_sha1']
            aa_length = transl_seq_meta[transl_key]['length']

            isoform_key = (gene_id, transc_id, transl_id, aa_length)
            isoform_meta[isoform_key] = {
                'gene_id': gene_id,
                'gene_version': gene_version,
                'transcript_id': transc_id,
                'transcript_version': transc_version,
                'translation_id': transl_id,
                'translation_version': transl_version,
                'translation_seq_sha1': transl_seq_sha1,
                'length': aa_length
            }

    return isoform_meta


def load_gencode_metadata(annot_file, transl_file):

    transl_seq_meta = SlottedDict()
    with open_as_text(transl_file) as f:
        for header, sequence in SimpleFastaParser(f):
            match = gencode_transl_hdr_regex.match(header)

            if is_par_y(match['ens_gene_id']):
                continue

            transl_id_ver = match['ens_transl_id']
            transl_seq_meta[transl_id_ver] = {
                'translation_seq_sha1': calc_sha1(sequence),
                'length': len(sequence)
            }

    transl_id_vers = set(transl_seq_meta.keys())

    annot_df = read_gtf(annot_file)

    annot_df = annot_df.loc[(annot_df['feature'] == 'transcript') &
                            annot_df['protein_id'].isin(transl_id_vers) &
                            ~annot_df['gene_id'].apply(is_par_y)]

    gene_to_isoforms = defaultdict(set)
    for gene_id_ver, transc_id_ver, transl_id_ver in zip(annot_df['gene_id'],
                                                         annot_df['transcript_id'],
                                                         annot_df['protein_id']):
        gene_to_isoforms[gene_id_ver].add((transc_id_ver, transl_id_ver))

    isoform_meta = SlottedDict()
    for gene_id_ver, isoforms in gene_to_isoforms.items():
        match = ens_id_regex.match(gene_id_ver)
        gene_id = match['bare_id']
        gene_version = match['id_version']

        for transc_id_ver, transl_id_ver in isoforms:
            match = ens_id_regex.match(transc_id_ver)
            transc_id = match['bare_id']
            transc_version = match['id_version']

            match = ens_id_regex.match(transl_id_ver)
            transl_id = match['bare_id']
            transl_version = match['id_version']

            transl_seq_sha1 = transl_seq_meta[transl_id_ver]['translation_seq_sha1']
            aa_length = transl_seq_meta[transl_id_ver]['length']

            isoform_key = (gene_id, transc_id, transl_id, aa_length)
            isoform_meta[isoform_key] = {
                'gene_id': gene_id,
                'gene_version': gene_version,
                'transcript_id': transc_id,
                'transcript_version': transc_version,
                'translation_id': transl_id,
                'translation_version': transl_version,
                'translation_seq_sha1': transl_seq_sha1,
                'length': aa_length
            }

    return isoform_meta


@contextmanager
def open_as_text(file, mode='r', **kwargs):
    if mode not in ('a', 'r', 'w', 'x'):
        raise ValueError(f"unknown file mode: '{mode}'")
    text_mode = mode + 't'
    if str(file).lower().endswith('.gz'):
        opener = functools.partial(gzip.open, mode=text_mode, **kwargs)
    else:
        opener = functools.partial(open, mode=text_mode, **kwargs)
    file_obj = None
    try:
        file_obj = opener(file)
        yield file_obj
    finally:
        if file_obj is not None:
            file_obj.close()


def resolve_source_slug(source_name, source_version):
    source_key = source_name.lower()
    try:
        source_prefix = src_name_to_prefix[source_key]
    except KeyError:
        raise ValueError(f"unknown source name: '{source_name}'")

    # We accept slight variations in GENCODE version numbers,
    # but we resolve these to the format expected by TRIFID.
    if source_key == 'gencode':
        match = gencode_version_regex.match(source_version)
        source_version = match['version']

    return f'{source_prefix}{source_version}'


ens_id_regex = re.compile('^(?P<bare_id>.+)\.(?P<id_version>[0-9]+)$')

exp_trifid_cols = set([
    'gene_id',
    'gene_name',
    'transcript_id',
    'translation_id',
    'flags',
    'ccdsid',
    'appris',
    'ann_type',
    'length',
    'trifid_score',
    'norm_trifid_score'
])

ensembl_transl_hdr_regex = re.compile(
    '^(?P<ens_transl_id>\S+)'
    ' pep (?P<location>\S+)'
    ' gene:(?P<ens_gene_id>\S+)'
    ' transcript:(?P<ens_transc_id>\S+)'
    ' gene_biotype:(?P<gene_biotype>\S+)'
    ' transcript_biotype:(?P<transc_biotype>\S+)'
    '(?: gene_symbol:(?P<gene_name>\S+))?'
    '(?: description:(?P<desc>.+))?'
)

gencode_transl_hdr_regex = re.compile(
    '^(?P<ens_transl_id>[^\|]+)\|'
    '(?P<ens_transc_id>[^\|]+)\|'
    '(?P<ens_gene_id>[^\|]+)\|'
    '(?P<hav_gene_id>[^\|]+)\|'
    '(?P<hav_transc_id>[^\|]+)\|'
    '(?P<transc_name>[^\|]+)\|'
    '(?P<gene_name>[^\|]+)\|'
    '(?P<aa_length>\d+)'
)

gencode_version_regex = re.compile('^v?M?(?P<version>[0-9]+)$')

load_metadata = OrderedDict([
    ('ensembl', load_ensembl_metadata),
    ('gencode', load_gencode_metadata)
])

out_trifid_cols = [
    'gene_id',
    'gene_version',
    'gene_name',
    'transcript_index',
    'transcript_id',
    'transcript_version',
    'translation_id',
    'translation_version',
    'translation_seq_sha1',
    'flags',
    'ccdsid',
    'appris',
    'ann_type',
    'length',
    'trifid_score',
    'norm_trifid_score'
]

src_name_to_prefix = OrderedDict([
    ('ensembl', 'e'),
    ('gencode', 'g')
])

supported_sources = list(src_name_to_prefix.keys())


TrifidRecord = namedtuple('TrifidRecord', out_trifid_cols)


ap = ArgumentParser(description=__doc__)
ap.add_argument('--assembly', required=True,
                help='Genome assembly used as input to TRIFID.')
ap.add_argument('--source-name', required=True,
                choices=supported_sources,
                help='Name of data source used as input to TRIFID.')
ap.add_argument('--source-version', required=True,
                help='Version of data source used as input to TRIFID.')
ap.add_argument('--input-file', metavar='PATH',
                help='If given, read TRIFID prediction'
                     ' data from the specified file.')
ap.add_argument('--annot-file', metavar='PATH', required=True,
                help='Input feature annotation file.')
ap.add_argument('--transl-file', metavar='PATH', required=True,
                help='Input feature translation sequence file.')
ap.add_argument('--output-dir', metavar='PATH', default=os.getcwd(),
                help='Output directory in which to save the'
                     ' processed TRIFID prediction file.')
args = ap.parse_args()

assembly = args.assembly
source_name = args.source_name
source_version = args.source_version
input_file = args.input_file
annot_file = args.annot_file
transl_file = args.transl_file
output_dir = args.output_dir


print('Resolving source slug …')
source_slug = resolve_source_slug(source_name, source_version)

print('Loading isoform metadata …')
isoform_meta = load_metadata[source_name](annot_file, transl_file)

with TemporaryDirectory() as tmp_dir:

    pred_file_path = os.path.join(tmp_dir, 'trifid_predictions.tsv.gz')
    if input_file is not None:
        print('Copying local TRIFID file …')
        pred_file_version = copy_local_trifid_file(input_file, pred_file_path)
    else:
        print('Downloading TRIFID file from GitLab …')
        pred_file_version = download_gitlab_trifid_file(assembly, source_slug,
                                                        pred_file_path)

    text_file_name = f'trifid_{assembly}_{source_slug}_{pred_file_version}.tsv'
    out_file_name = f'{text_file_name}.gz'

    print('Reading TRIFID predictions …')
    gene_to_recs = defaultdict(list)
    with open_as_text(pred_file_path) as f:
        reader = csv.DictReader(f, delimiter='\t', lineterminator='\n')

        obs_trifid_cols = set(reader.fieldnames)
        missing_cols = exp_trifid_cols - obs_trifid_cols
        if missing_cols:
            raise ValueError(f"missing column(s): {', '.join(missing_cols)}")

        unknown_cols = obs_trifid_cols - exp_trifid_cols
        if unknown_cols:
            raise ValueError(f"unknown column(s): {', '.join(unknown_cols)}")

        for rec in reader:
            gene_id = rec['gene_id']
            transc_id = rec['transcript_id']
            transl_id = rec['translation_id']

            try:
                aa_length = int(rec['length'])
            except ValueError as e:
                # Accept a float if it is an integer.
                aa_length = float(rec['length'])
                if aa_length.is_integer():
                    aa_length = round(aa_length)
                else:
                    raise e

            transc_key = (gene_id, transc_id, transl_id, aa_length)
            meta_rec = isoform_meta[transc_key]

            for k in ('gene_id', 'gene_version',
                      'transcript_id', 'transcript_version',
                      'translation_id', 'translation_version',
                      'translation_seq_sha1', 'length'):
                rec[k] = meta_rec[k]

            gene_to_recs[gene_id].append(rec)

    out_file_cols = out_trifid_cols.copy()
    s_idx = 1 + out_file_cols.index('gene_id')
    p_idx = 1 + out_file_cols.index('transcript_index')

    print('Sorting TRIFID predictions …')
    text_file_path = os.path.join(tmp_dir, text_file_name)
    with open_as_text(text_file_path, mode='w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')

        # Prepend hash symbol so tabix will take first line as header.
        out_file_cols[0] = f'#{out_file_cols[0]}'
        writer.writerow(out_file_cols)

        for gene_id in sorted(gene_to_recs.keys()):
            recs = gene_to_recs[gene_id]
            recs.sort(
                key=lambda x: (-Decimal(x['trifid_score']), x['transcript_id'])
            )
            for transc_idx, rec in enumerate(recs, start=1):
                rec['transcript_index'] = transc_idx

            writer.writerows([TrifidRecord(**rec) for rec in recs])

    print('Compressing TRIFID file …')
    bgzip_file_path = os.path.join(tmp_dir, out_file_name)
    bgzip_cmd_args = ['bgzip', text_file_path]
    subprocess.run(bgzip_cmd_args, check=True)

    print('Indexing TRIFID file …')
    out_index_name = f'{out_file_name}.tbi'
    tmp_index_path = os.path.join(tmp_dir, out_index_name)
    out_index_path = os.path.join(output_dir, out_index_name)
    tabix_cmd_args = [
        'tabix',
        '-s', str(s_idx),
        '-b', str(p_idx),
        '-e', str(p_idx),
        bgzip_file_path
    ]
    subprocess.run(tabix_cmd_args, check=True)

    out_file_path = os.path.join(output_dir, out_file_name)
    print(f"Saving prepared TRIFID file: '{out_file_path}' …")
    os.makedirs(output_dir, exist_ok=True)
    shutil.move(bgzip_file_path, out_file_path)
    shutil.move(tmp_index_path, out_index_path)

print('Done.')
