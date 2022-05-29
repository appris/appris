#!/usr/bin/env python3
"""Prepare PDBsum BLAST database for Matador3D."""

from argparse import ArgumentParser
from datetime import datetime, timezone
import os
import re
import shutil
from subprocess import PIPE, run
from tempfile import TemporaryDirectory

from Bio.SeqIO.FastaIO import SimpleFastaParser
import requests


ap = ArgumentParser(description=__doc__)
ap.add_argument('--output-dir', metavar='PATH', default=os.getcwd(),
                help='Output directory in which to save'
                     ' the PDBsum BLAST database files.')

args = ap.parse_args()
output_dir = args.output_dir


curr_pdb_id_url = 'https://data.rcsb.org/rest/v1/holdings/current/entry_ids'
in_file_url = 'https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/pdblib.fasta'

print('Fetching current PDB ID set...')
with requests.get(curr_pdb_id_url) as r:
    r.raise_for_status()
    curr_pdb_ids = set(r.json())

with TemporaryDirectory() as tmp_dir:

    print('Downloading PDBsum FASTA file...')
    tmp_fasta_file = os.path.join(tmp_dir, 'pdblib.fasta')
    with requests.get(in_file_url, stream=True) as r:
        with open(tmp_fasta_file, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    datestamp = datetime.utcnow().strftime('%Y%m%d')
    blastdb_name = f'pdbsum_{datestamp}'
    blastdb_path = os.path.join(output_dir, blastdb_name)
    out_fasta_file = os.path.join(output_dir, f'{blastdb_name}.fasta')

    print('Filtering PDBsum FASTA file...')
    pdbsum_seqid_regex = re.compile('^(?P<pdb_id>[1-9][a-z0-9]{3}):'
                                    '(?P<chain_id>[A-Za-z0-9]{1,4})?$')
    num_entries = 0
    obsolete = set()
    chainless = set()
    with open(tmp_fasta_file, 'r') as in_f, open(out_fasta_file, 'w') as out_f:
        for header, sequence in SimpleFastaParser(in_f):
            num_entries += 1

            seqid, *other_text = header.split(maxsplit=1)
            match = pdbsum_seqid_regex.match(seqid)

            # We capitalise the bare PDB ID to match the
            # format of the current IDs returned by PDB.
            pdb_id = match['pdb_id'].upper()

            # Chain IDs are case-sensitive, so we preserve case.
            chain_id = match['chain_id']

            dropping = False
            if pdb_id not in curr_pdb_ids:
                obsolete.add(seqid)
                dropping = True
            if chain_id is None:
                chainless.add(seqid)
                dropping = True

            if dropping:
                continue

            seqid = f'PDB:{pdb_id}_{chain_id}'
            header = ' '.join([seqid] + other_text)
            out_f.write(f'>{header}\n{sequence}\n')

    num_dropped = len(obsolete | chainless)
    print(f'Processed {num_entries} PDBsum entries')
    print(f'Found {len(obsolete)} obsolete entries...')
    print(f'Found {len(chainless)} entries lacking a chain ID...')
    print(f'Dropped {num_dropped} entries...')

    print(f"Making BLAST database '{blastdb_path}'...")
    blastdb_log = os.path.join(tmp_dir, 'formatdb.log')
    cmd_args = [
        'formatdb',
        '-i', out_fasta_file,
        '-p', 'T',
        '-l', blastdb_log,
        '-t', blastdb_name,
        '-n', blastdb_path
    ]
    run(cmd_args,stdout=PIPE,stderr=PIPE,check=True)

print('Done.')
