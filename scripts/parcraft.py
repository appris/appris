#!/usr/bin/env python3
"""Utilities for handling pseudoautosomal (PAR) features."""

from argparse import ArgumentParser
from collections import defaultdict, namedtuple, OrderedDict
from contextlib import closing, contextmanager
import functools
import gzip
import logging
import os
import re
import sys
from tempfile import TemporaryDirectory

import gffutils
from pysam import FastaFile


__all__ = ['filter_par_features', 'tag_par_features']


log_level_map = OrderedDict([
    ('debug', logging.DEBUG),
    ('info', logging.INFO),
    ('warning', logging.WARNING),
    ('error', logging.ERROR),
    ('critical', logging.CRITICAL)
])

logger = logging.Logger(__name__)


Feature = namedtuple('Feature', ['seqid', 'source', 'featuretype', 'start', 'end',
                                   'score', 'strand', 'frame', 'attributes'])


def filter_par_features(in_data_file, out_data_file):
    """Filter PAR features from a GTF/GFF3 file."""

    logger.info('Filtering PAR features…')
    with open_as_text(in_data_file, 'r') as in_f:
        with open_as_text(out_data_file, 'w') as out_f:
            for line in in_f:
                if not line.startswith('#'):
                    feature, _ = parse_dataline(line)
                    if ('tag' in feature.attributes and
                            'PAR' in feature.attributes['tag']):
                        continue
                out_f.write(line)

    logger.info('Done.')


def genbank_genes(db, *args, **kwargs):
    for feature in db.all_features(*args, **kwargs):
        if ('gbkey' in feature.attributes and
                feature.attributes['gbkey'][0] == 'Gene'):
            yield feature


def get_ncbi_gene_id(feature):
    dbxrefs = parse_dbxrefs(feature)
    try:
        return dbxrefs['GeneID']
    except KeyError:
        return None


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


def parse_dataline(line):
    """Parse data line to a feature."""

    line = line.rstrip('\r\n')
    fields = line.split('\t')
    try:
        (chrom, source, feat_type, start, end,
         score, strand, frame, attr_field) = fields
    except ValueError:
        raise ValueError(f'expected 9 tab-delimited fields, not {len(fields)}')

    attrs = OrderedDict()
    inferred_formats = defaultdict(int)
    if attr_field != '.':
        if attr_field.endswith(';'):
            attr_field = attr_field[:-len(';')]
        for subfield in re.split('\s*;\s*', attr_field):
            match = re.match('^\s*(?P<key>\S+) "?(?P<value>[^"]*)"?\s*$', subfield)
            if match:
                inferred_formats['gtf'] += 1
            else:
                match = re.match('^\s*(?P<key>[^=]+)=(?P<value>.*)\s*$', subfield)
                if match:
                    inferred_formats['gff3'] += 1
                else:
                    raise ValueError(f"cannot parse attribute: '{subfield}'")
            key = match['key']
            values = match['value'].split(',')
            try:
                attrs[key].extend(values)
            except KeyError:
                attrs[key] = values

    feature = Feature(chrom, source, feat_type, start,
                      end, score, strand, frame, attrs)

    feat_fmt = None
    if inferred_formats:
        feat_fmt = max((fmt for fmt in inferred_formats.keys()),
                       key=lambda k: inferred_formats[k])

    return feature, feat_fmt


def parse_dbxrefs(feature):
    return dict([x.split(':', maxsplit=1)
                 for x in feature.attributes.get('Dbxref', [])])


def serialise_feature(feature, feature_format):
    """Serialise feature to a data line."""
    fields = list(feature)
    attrs = fields.pop()

    if feature_format == 'gtf':
        attr_field = ';'.join(f'''{key} "{','.join(values)}"'''
                              for key, values in attrs.items())
    elif feature_format == 'gff3':
        attr_field = ';'.join(f"{key}={','.join(values)}"
                              for key, values in attrs.items())
    else:
        raise ValueError(f"unknown feature format: '{feature_format}'")

    return '{}\n'.format('\t'.join(fields + [attr_field]))


def tag_par_features(in_data_file, asm_conf_file, out_data_file, transl_file=None):
    """
    Add PAR tag to PAR features in a RefSeq GFF3 file.

    Groups of genes are assigned a PAR tag if they are located in
    a known PAR region, or if every gene in the group is located
    on an allosome and all genes in the group have the same set
    of translation sequences.

    Note that the aim is to identify and tag apparently duplicated
    PAR genes, and not to identify PAR genes as such.
    """

    logger.info('Loading assembly metadata…')
    asm_db = gffutils.create_db(asm_conf_file, ':memory:')

    with TemporaryDirectory() as tmp_dir:

        logger.info('Loading input data…')
        data_db_file = os.path.join(tmp_dir, 'gff.db')
        data_db = gffutils.create_db(in_data_file, data_db_file,
                                     merge_strategy='create_unique')

        if transl_file is not None:
            tmp_transl_file = os.path.join(tmp_dir, 'transl.fa')
            transl_file_path = os.path.abspath(transl_file)
            os.symlink(transl_file_path, tmp_transl_file)

        logger.info('Getting metadata on chromosomes…')
        chrom_id_map = dict()
        chrom_regions = list()
        known_chroms = set()
        allosomes = set()
        for feature in asm_db.features_of_type('chromosome'):
            dbxrefs = parse_dbxrefs(feature)
            try:
                refseq_acc = dbxrefs['RefSeq']
            except KeyError:
                raise ValueError(
                    f"RefSeq accession not found for chromosome '{feature.seqid}'")
            chrom_id_map[refseq_acc] = feature.attributes['Name'][0]

            chrom_region = (refseq_acc, feature.start, feature.end)
            chrom_regions.append(chrom_region)
            known_chroms.add(refseq_acc)

            tags = feature.attributes.get('tag', [])
            if 'allosome' in tags:
                allosomes.add(refseq_acc)

        logger.info('Getting metadata on known PAR regions…')
        known_par_regions = list()
        for feature in asm_db.features_of_type('region'):
            tags = feature.attributes.get('tag', [])
            if 'PAR' in tags:
                parent_id = feature.attributes['Parent'][0]
                parent = asm_db[parent_id]
                dbxrefs = parse_dbxrefs(parent)
                try:
                    refseq_acc = dbxrefs['RefSeq']
                except KeyError:
                    raise ValueError(
                        f"RefSeq accession not found for chromosome '{parent.seqid}'")
                region = (refseq_acc, feature.start, feature.end)
                known_par_regions.append(region)

        logger.info('Searching for duplicate NCBI gene IDs…')
        gene_to_chroms = defaultdict(set)
        for chrom_region in chrom_regions:
            chrom, start, end = chrom_region
            for gene in genbank_genes(data_db, limit=chrom_region):
                gene_id = get_ncbi_gene_id(gene)
                if gene_id is None:
                    raise ValueError(f"cannot find GeneID in gene feature: '{gene}'")
                gene_to_chroms[gene_id].add(chrom)

        dup_gene_ids = set([gene_id for gene_id, chroms in gene_to_chroms.items()
                            if len(chroms) > 1])

        # In flagging duplicates, we look for groups of genes having the same
        # NCBI GeneID and located on different autosomes. We do not treat
        # as duplicates pairs of genes located on allosomes, or those
        # groups of genes with the same NCBI GeneID on the same chromosome…
        # "…because they are considered to be different parts or alleles of the same gene."
        # — https://www.ncbi.nlm.nih.gov/datasets/docs/about-ncbi-gff3/
        cand_par_gene_ids = set()
        for gene_id in dup_gene_ids:
            gene_chroms = gene_to_chroms[gene_id]
            if gene_chroms <= allosomes:
                cand_par_gene_ids.add(gene_id)
            else:
                gene_chrom_text = "'{}'".format("', '".join(gene_chroms))
                logger.warning(f"duplicate gene ID '{gene_id}' found"
                               f" on chromosomes: {gene_chrom_text}")

        logger.info('Assembling set of genes in known PAR regions')
        known_par_gene_ids = set()
        for par_region in known_par_regions:
            for gene in genbank_genes(data_db, limit=par_region):
                gene_id = get_ncbi_gene_id(gene)
                if gene_id in cand_par_gene_ids:
                    known_par_gene_ids.add(gene_id)

        cand_par_gene_ids -= known_par_gene_ids

        inferred_par_gene_ids = set()
        if cand_par_gene_ids and transl_file is not None:
            logger.info('Inferring PAR genes from duplicates in allosomes…')

            cand_par_gene_features = defaultdict(list)
            for chrom_region in chrom_regions:
                for gene in genbank_genes(data_db, limit=chrom_region):
                    gene_id = get_ncbi_gene_id(gene)
                    if gene_id in cand_par_gene_ids:
                        cand_par_gene_features[gene_id].append(gene)
            cand_par_gene_features = dict(cand_par_gene_features)

            for gene_id, genes in cand_par_gene_features.items():

                transl_seq_sets = set()
                for gene in genes:

                    transl_ids = set([
                        child.attributes['protein_id'][0]
                        for child in data_db.children(gene)
                        if 'protein_id' in child.attributes
                    ])

                    transl_seq_set = set()
                    with closing(FastaFile(tmp_transl_file)) as f:
                        for transl_id in transl_ids:
                            transl_seq = f.fetch(transl_id)
                            transl_seq_set.add(transl_seq)

                    transl_seq_sets.add(frozenset(transl_seq_set))

                if len(transl_seq_sets) == 1:
                    for gene in genes:
                        gene_id = get_ncbi_gene_id(gene)
                        inferred_par_gene_ids.add(gene_id)

        par_gene_ids = known_par_gene_ids | inferred_par_gene_ids

    logger.info('Tagging PAR features…')
    with open_as_text(in_data_file, 'r') as in_f:
        with open_as_text(out_data_file, 'w') as out_f:
            for line in in_f:
                if not line.startswith('#'):
                    feature, feat_fmt = parse_dataline(line)
                    canonic_chrom_id = chrom_id_map.get(feature.seqid)
                    if canonic_chrom_id == 'Y':  # TODO: eliminate hard-coded value
                        gene_id = get_ncbi_gene_id(feature)
                        if gene_id in par_gene_ids and not (
                                'tag' in feature.attributes and
                                'PAR' in feature.attributes['tag']):
                            if gene_id not in known_par_gene_ids:
                                logger.info(
                                    f"tagging feature of inferred PAR gene '{gene_id}'")
                            tags = feature.attributes.setdefault('tag', [])
                            tags.append('PAR')
                            line = serialise_feature(feature, feature_format=feat_fmt)
                out_f.write(line)

    logger.info('Done.')


if __name__ == '__main__':

    parser = ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(title='commands', dest='command')

    filter_parser = subparsers.add_parser('filter',
                                          description=filter_par_features.__doc__,
                                          help='filter PAR features')
    filter_parser.add_argument('--in-data-file', required=True,
                               help='Input annotation GTF/GFF3 file.')
    filter_parser.add_argument('--out-data-file', required=True,
                               help='Output annotation file with PAR features'
                                    ' removed from the Y chromosome.')
    filter_parser.set_defaults(func=filter_par_features)

    tag_parser = subparsers.add_parser('tag', description=tag_par_features.__doc__,
                                       help='tag PAR features')
    tag_parser.add_argument('--in-data-file', required=True,
                            help='Input RefSeq GFF3 file.')
    tag_parser.add_argument('--asm-conf-file', required=True,
                            help='Assembly metadata GFF3 file.')
    tag_parser.add_argument('--out-data-file', required=True,
                            help='Output RefSeq GFF3 annotation file,'
                                 ' with PAR tags added as appropriate.')
    tag_parser.add_argument('--transl-file',
                            help='Optional input translation sequence file, which is used to'
                                 ' infer PAR genes in addition to those in known PAR regions.')
    tag_parser.set_defaults(func=tag_par_features)

    parser.add_argument('--log-level', choices=log_level_map.keys(),
                        help='log level')

    args = parser.parse_args()
    kwargs = vars(args)
    log_level = kwargs.pop('log_level')
    command = kwargs.pop('command')
    func = kwargs.pop('func', None)

    handler = logging.StreamHandler(sys.stderr)
    if log_level is not None:
        logger.setLevel(log_level_map[log_level])
        handler.setLevel(log_level_map[log_level])
    formatter = logging.Formatter('[%(asctime)s - %(filename)s] %(message)s',
                                  datefmt='%Y-%m-%dT%H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if command is not None:
        func(**kwargs)
    else:
        parser.print_help()
