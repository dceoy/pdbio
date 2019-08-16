#!/usr/bin/env python

import logging

from .vcfdataframe import VcfDataFrame


def identify_variants(fa_path, chrom, variants):
    """
    Args:
        fa_path (str): path to a reference genome FASTA file
        chrom (str): chromosome name
        variants (list): variant expressions
    """
    logger = logging.getLogger(__name__)
    assert len(set(variants)) == len(variants), 'Duplicated arguments'
    vars = [_variant2dict(variant=s, chrom=chrom) for s in variants]
    logger.debug('vars: {}'.format(vars))
    edges = [min([v['pos'] for v in vars]), max([v['endpos'] for v in vars])]
    logger.debug('edges: {}'.format(edges))
    chrom_seq = _read_chrom_seq(fa_path=fa_path, chrom=chrom).upper()
    ref_seq = chrom_seq[(edges[0] - 1):edges[1]]
    print('> {0}:{1:d}-{2:d}'.format(chrom, *edges), flush=True)
    ref_id = 'FASTA sequence'
    stdout_str = (
        '  {0:<' + str(max(len(ref_id), *[len(v['id']) for v in vars]))
        + '} : \t{1}'
    )
    logger.debug('stdout_str: {}'.format(stdout_str))
    print(stdout_str.format(ref_id, ref_seq), flush=True)
    id_seq_set = set()
    for v in vars:
        logger.debug('v: {}'.format(v))
        e_seqs = [
            chrom_seq[(edges[0] - 1):(v['pos'] - 1)],
            chrom_seq[v['endpos']:edges[1]]
        ]
        logger.debug('e_seqs: {}'.format(e_seqs))
        bases = [(e_seqs[0] + s + e_seqs[1]) for s in [v['ref'], v['alt']]]
        logger.info('bases (REF -> ALT): {0} -> {1}'.format(*bases))
        assert bases[0] == ref_seq, 'REF bases discord from the genome'
        id_seq_set.add(bases[1])
        print(stdout_str.format(v['id'], bases[1]), flush=True)
    print('Detected variants :\t{}'.format(len(id_seq_set)))


def _variant2dict(variant, chrom):
    logger = logging.getLogger(__name__)
    pos_str, ref_alt = tuple(variant.split(':', maxsplit=1))
    logger.debug('{0} -> {1}, {2}'.format(variant, pos_str, ref_alt))
    pos = int(pos_str)
    ref, alt = tuple(ref_alt.split('_', maxsplit=1))
    endpos = pos + len(ref) - 1
    return {
        'id': '{0}:{1:d}-{2:d}:{3}>{4}'.format(chrom, pos, endpos, ref, alt),
        'pos': pos, 'endpos': endpos, 'ref': ref.upper(), 'alt': alt.upper()
    }


def _read_chrom_seq(fa_path, chrom):
    logger = logging.getLogger(__name__)
    logger.info('Read {0} in FASTA: {1}'.format(chrom, fa_path))
    vdf = VcfDataFrame()
    chrom_seq = ''
    with vdf.open_readable_file(path=vdf.normalize_path(path=fa_path)) as f:
        loading = False
        for s in f:
            if s.startswith('>'):
                fa_chrom = s.strip()[1:]
                loading = (fa_chrom == chrom)
                if loading or not chrom_seq:
                    logger.debug(('Load ' if loading else 'Skip ') + fa_chrom)
                else:
                    logger.debug('len(chrom_seq): {}'.format(len(chrom_seq)))
                    break
            elif loading:
                chrom_seq += s.strip()
    assert bool(chrom_seq), '{0} not found: {1}'.format(chrom, fa_path)
    return chrom_seq
