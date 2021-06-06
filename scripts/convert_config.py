#!/usr/bin/env python3
'''Convert APPRIS config file from Excel to JSON.

Because (some) people prefer Excel,
while machines prefer JSON.
'''

from argparse import ArgumentParser
from collections import OrderedDict
import json

from openpyxl import load_workbook
import pandas as pd


def comma_separated_list(s):
    value = str(s).strip()
    if s:
        value = tuple(x.strip() for x in value.split(','))
    else:
        value = None
    return value

def nullable_bool(s):
    value = str(s).strip()
    if value == 'true':
        value = True
    elif value == 'false':
        value = False
    elif value == '':
        value = None
    else:
        raise ValueError(f"cannot infer bool value from string: '{s}'")
    return value

def nullable_integer(s):
    value = str(s).strip()
    if value != '':
        value = int(value)
    else:
        value = None
    return value

def nullable_string(s):
    value = str(s).strip()
    if value == '':
        value = None
    return value


converters = {
    'species': {
        'id': str,
        'order': int,
        'scientific': str,
        'taxid': nullable_integer,
        'common': str,
        'group': str,
        'category': nullable_string,
        'official': str
    },
    'assemblies': {
        'species': str,
        'id': str,
        'name': str
    },
    'datasets': {
        'assembly': str,
        'id': str,
        'type': str,
        'queryable': nullable_bool,
        'source_name': str,
        'source_label': str,
        'source_version': str,
        'database_name': nullable_string,
        'database_inifile': nullable_string,
        'datafiles_principal_isoforms': comma_separated_list,
        'datafiles_appris_scores': comma_separated_list,
        'datafiles_functional_residues': comma_separated_list,
        'datafiles_tertiary_structure': comma_separated_list,
        'datafiles_tertiary_structure_2': comma_separated_list,
        'datafiles_vertebrates_conservation': comma_separated_list,
        'datafiles_whole_domains': comma_separated_list,
        'datafiles_transmembrane_helices': comma_separated_list,
        'datafiles_peptides': comma_separated_list,
        'pipeline_envfile': nullable_string,
        'trifid_release': nullable_string
    },
    'datafiles': {
        'id': str,
        'label': str,
        'desc': str,
        'file': nullable_string
    }
}

dset_key_maps = OrderedDict([
    ('source', OrderedDict([
        ('name', 'source_name'),
        ('label', 'source_label'),
        ('version', 'source_version')
    ])),
    ('database', OrderedDict([
        ('name', 'database_name'),
        ('inifile', 'database_inifile'),
    ])),
    ('datafiles', OrderedDict([
        ('principal_isoforms', 'datafiles_principal_isoforms'),
        ('appris_scores', 'datafiles_appris_scores'),
        ('functional_residues', 'datafiles_functional_residues'),
        ('tertiary_structure', 'datafiles_tertiary_structure'),
        ('tertiary_structure_2', 'datafiles_tertiary_structure_2'),
        ('vertebrates_conservation', 'datafiles_vertebrates_conservation'),
        ('whole_domains', 'datafiles_whole_domains'),
        ('transmembrane_helices', 'datafiles_transmembrane_helices'),
        ('peptides', 'datafiles_peptides')
    ])),
    ('pipeline', OrderedDict([
        ('envfile', 'pipeline_envfile')
    ])),
    ('trifid', OrderedDict([
        ('release', 'trifid_release')
    ]))
])


ap = ArgumentParser(description=__doc__)
ap.add_argument('-i', '--input-file', metavar='PATH', required=True,
                help='APPRIS config Excel file')
ap.add_argument('-o', '--output-file', metavar='PATH', required=True,
                help='APPRIS config JSON file')
args = ap.parse_args()

input_file = args.input_file
output_file = args.output_file


config = OrderedDict()

wb = load_workbook(input_file)

version_df = pd.read_excel(input_file, 'version',
                           dtype=object, na_filter=False,
                           index_col=False, engine='openpyxl')
config['version'] = version_df.loc[0, 'version']

df = dict()
for sheet_name in ('species', 'assemblies', 'datasets', 'datafiles'):
    df[sheet_name] = pd.read_excel(input_file, sheet_name,
        converters=converters[sheet_name], na_filter=False,
        index_col=False, engine='openpyxl')

asm_to_specie = dict()
specie_to_asms = OrderedDict()
for specie_id, asm_id in zip(df['assemblies'].species, df['assemblies'].id):
    try:
        specie_to_asms[specie_id].append(asm_id)
    except KeyError:
        specie_to_asms[specie_id] = [asm_id]
    asm_to_specie[asm_id] = specie_id

dset_uids = list()
asm_to_dsets = OrderedDict()
for asm_id, dset_id in zip(df['datasets'].assembly, df['datasets'].id):
    # Dataset IDs are reused in different species, so we
    # create a dataset UID by prepending the species ID.
    dset_uid = f'{asm_to_specie[asm_id]}_{dset_id}'
    try:
        asm_to_dsets[asm_id].append(dset_uid)
    except KeyError:
        asm_to_dsets[asm_id] = [dset_uid]
    dset_uids.append(dset_uid)

# Set DataFrame index to dataset UIDs.
df['datasets'].index = dset_uids
if df['datasets'].duplicated().any():
    raise ValueError('dataset IDs must be unique within a species')

for sheet_name in ('species', 'assemblies', 'datafiles'):
    df[sheet_name].set_index('id', drop=False, inplace=True,
        verify_integrity=True)

config['species'] = OrderedDict()
for _, specie_row in df['species'].iterrows():
    specie_id = specie_row['id']

    specie_cfg = OrderedDict()
    specie_cfg['order'] = specie_row['order']
    specie_cfg['scientific'] = specie_row['scientific']
    if 'taxid' in specie_row and specie_row['taxid'] is not None:
        specie_cfg['taxid'] = specie_row['taxid']
    specie_cfg['common'] = specie_row['common']
    specie_cfg['group'] = specie_row['group']
    if specie_row['category'] is not None:
        specie_cfg['category'] = specie_row['category']
    specie_cfg['official'] = specie_row['official']

    specie_cfg['assemblies'] = list()
    for asm_id in specie_to_asms[specie_id]:
        asm_row = df['assemblies'].loc[asm_id, :]

        asm_cfg = OrderedDict()
        asm_cfg['id'] = asm_row['id']
        asm_cfg['name'] = asm_row['name']
        asm_cfg['datasets'] = list()

        for dset_uid in asm_to_dsets[asm_id]:
            dset_row = df['datasets'].loc[dset_uid, :]

            dset_cfg = OrderedDict()
            dset_cfg['id'] = dset_row['id']
            dset_cfg['type'] = dset_row['type']

            if ('queryable' in dset_row and
                    dset_row['queryable'] is not None):
                dset_cfg['queryable'] = bool(dset_row['queryable'])

            for section, key_pairs in dset_key_maps.items():
                rel_key_pairs = [(k, x) for k, x in key_pairs.items()
                                 if dset_row[x] is not None]
                if rel_key_pairs:
                    dset_cfg[section] = OrderedDict([
                        (key, dset_row[col_name])
                        for key, col_name in rel_key_pairs
                    ])

            asm_cfg['datasets'].append(dset_cfg)
        specie_cfg['assemblies'].append(asm_cfg)
    config['species'][specie_id] = specie_cfg

config['datafiles'] = list()
for _, row in df['datafiles'].iterrows():
    config['datafiles'].append(OrderedDict([
        (k, row[k]) for k in ('id', 'label', 'desc', 'file')
        if row[k] is not None
    ]))

with open(output_file, 'w') as f:
    json.dump(config, f, indent=4)
