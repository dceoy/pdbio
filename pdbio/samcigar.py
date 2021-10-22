#!/usr/bin/env python
"""
SAM CIGAR Parser.
https://github.com/dceoy/pdbio
"""

import re

import numpy as np
import pandas as pd


def cigar2reflen(cigar):
    """
    Args:
        cigar (str): CIGAR string of SAM

    Returns:
        int: length of a consumed reference bases

    Examples:
        >>> cigar2reflen(cigar='26S15M4D32M3I75M')
        126
    """
    return np.array([
        (s or '0')
        for s in re.split('M|D|N|=|X', re.sub(r'[0-9]+[ISHP]', '', cigar))
    ]).astype(np.int16).sum()


def cigar2qlen(cigar):
    """
    Args:
        cigar (str): CIGAR string of SAM

    Returns:
        int: length of a consumed query bases

    Examples:
        >>> cigar2qlen(cigar='26S15M4D32M3I75M')
        151
    """
    return np.array([
        (s or '0')
        for s in re.split('M|I|S|=|X', re.sub(r'[0-9]+[DNHP]', '', cigar))
    ]).astype(np.int16).sum()


def cigar2oplen(cigar):
    """
    Args:
        cigar (str): CIGAR string of SAM

    Returns:
        dict: total lengths of CIGAR operations

    Examples:
        >>> cigar2oplen(cigar='26S15M4D32M3I75M')
        {'D': 4, 'I': 3, 'M': 122, 'S': 26}
    """
    return pd.DataFrame(
        [
            (k, int(v)) for k, v in zip(
                re.sub(r'[0-9]', '', cigar),
                re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
            )
        ],
        columns=['op', 'len']
    ).groupby('op')['len'].sum().to_dict()


def cigar2chrs(cigar, only_aligned=False):
    """
    Args:
        cigar (str): CIGAR string of SAM
        only_aligned (bool): cut unmapped characters on edges

    Returns:
        str: CIGAR character sequence

    Examples:
        >>> cigar2chrs(cigar='6S5M4D12M3I5M')
        'SSSSSSMMMMMDDDDMMMMMMMMMMMMIIIMMMMM'
        >>> cigar2chrs(cigar='6S5M4D12M3I5M', only_aligned=True)
        'MMMMMDDDDMMMMMMMMMMMMIIIMMMMM'
    """
    chrs = ''.join([
        (k * int(v)) for k, v in zip(
            re.sub(r'[0-9]', '', cigar),
            re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
        )
    ])
    if only_aligned:
        return re.sub(r'^[ISHP]+', '', re.sub(r'[ISHP]+$', '', chrs))
    else:
        return chrs


def md2chrs(md):
    """
    Args:
        md (str): MD tag value of SAM

    Returns:
        str: MD character sequence
            (matched: '=', unmatched: 'X', deleted: 'D')

    Examples:
        >>> md2chrs(md='5C5^T16')
        '=====X=====D================'
        >>> md2chrs(md='16^AGTTTCA33')
        '================DDDDDDD================================='
        >>> md2chrs(md='7A0A30C9')
        '=======XX==============================X========='
    """
    not_match_md = [
        (('D' * len(s[1:])) if s.startswith('^') else ('X' if s else s))
        for s in re.split(r'[0-9]+', md)
    ]
    n_not_match_md = len(not_match_md)
    match_md = [
        ('=' * int(s)) for s in re.split(r'[\^A-Z]+', md) if s
    ]
    if n_not_match_md > len(match_md):
        match_md.append('')
    assert n_not_match_md == len(match_md)
    return ''.join([(u + m) for m, u in zip(match_md, not_match_md)])


def cigar2matchchrs(cigar, md):
    """
    Args:
        cigar (str): CIGAR string of SAM
        md (str): MD tag value of SAM

    Returns:
        str: CIGAR+MD character sequence
            (matched: '=', unmatched: 'X', deleted: 'D', inserted: 'I')

    Examples:
        >>> cigar2matchchrs(cigar='11M1D16M', md='5C5^T16')
        '=====X=====D================'
        >>> cigar2matchchrs(cigar='13M1I3M7D33M', md='16^AGTTTCA33')
        '=============I===DDDDDDD================================='
        >>> cigar2matchchrs(cigar='49M', md='7A0A30C9')
        '=======XX==============================X========='
    """
    md_chrs = md2chrs(md=md)
    md_i = 0
    aligned_cigar_list = list(cigar2chrs(cigar=cigar, only_aligned=True))
    for i, c in enumerate(aligned_cigar_list):
        if c == 'M':
            aligned_cigar_list[i] = md_chrs[md_i]
            md_i += 1
        elif c in 'D=X':
            md_i += 1
    return ''.join(aligned_cigar_list)


def seq2alignseq(seq, cigar):
    """
    Args:
        seq (str): SEQ string of SAM
        cigar (str): CIGAR string of SAM

    Returns:
        str: aligned query sequence

    Examples:
        >>> seq2alignseq(seq='TACAGCAGACGGGACCTTTTTGGTA', cigar='3S20M2S')
        'AGCAGACGGGACCTTTTTGG'
    """
    ops = [
        (k, int(v)) for k, v in zip(
            re.sub(r'[0-9]', '', cigar),
            re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
        )
    ]
    aligned_seq = seq
    for k, v in ops:
        if k in 'ISHP':
            aligned_seq = aligned_seq[v:]
        else:
            break
    for k, v in reversed(ops):
        if k in 'ISHP':
            aligned_seq = aligned_seq[:-v]
        else:
            break
    return aligned_seq


if __name__ == '__main__':
    import doctest
    doctest.testmod()
