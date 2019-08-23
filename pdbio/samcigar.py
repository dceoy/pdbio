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
        sum (bool): sum lengths by character

    Returns:
        dict: total lengths of CIGAR operations

    Examples:
        >>> cigar2oplen(cigar='26S15M4D32M3I75M')
        {'D': 4, 'I': 3, 'M': 122, 'S': 26}
    """
    return pd.DataFrame(
        [
            (k, np.int16(v)) for k, v in zip(
                re.sub(r'[0-9]', '', cigar),
                re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
            )
        ],
        columns=['op', 'len']
    ).groupby('op')['len'].sum().to_dict()


def cigar2alnrange(cigar):
    """
    Args:
        cigar (str): CIGAR string of SAM

    Returns:
        tuple: aligned sequence range of read

    Examples:
        >>> cigar2alnrange(cigar='26S15M4D32M3I75M')
        (26, 150)
    """
    ops = [
        (k, np.int16(v)) for k, v in zip(
            re.sub(r'[0-9]', '', cigar),
            re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
        )
    ]
    seq_range = [0, sum([v for k, v in ops if k in 'MIS=X']) - 1]
    for k, v in ops:
        if k in 'ISHP':
            seq_range[0] += v
        else:
            break
    for k, v in reversed(ops):
        if k in 'ISHP':
            seq_range[1] -= v
        else:
            break
    return tuple(seq_range)


def cigar2chars(cigar, only_aligned=False):
    """
    Args:
        cigar (str): CIGAR string of SAM
        only_aligned (bool): cut unmapped characters on edges

    Returns:
        str: CIGAR character sequence

    Examples:
        >>> cigar2chars(cigar='6S5M4D12M3I5M')
        'SSSSSSMMMMMDDDDMMMMMMMMMMMMIIIMMMMM'
    """
    chars = ''.join([
        (k * int(v)) for k, v in zip(
            re.sub(r'[0-9]', '', cigar),
            re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
        )
    ])
    if only_aligned:
        return re.sub(r'^[ISHP]+', '', re.sub(r'[ISHP]+$', '', chars))
    else:
        return chars


def cigar2matchchars(cigar, md):
    """
    Args:
        cigar (str): CIGAR string of SAM
        md (str): MD tag value of SAM

    Returns:
        str: CIGAR character sequence
            (match: '=', unmatch: 'X', deletion: 'D')

    Examples:
        >>> cigar2matchchars(cigar='11M1D16M', md='5C5^T16')
        '=====X=====D================'
    """
    aln_cigar_list = list(
        re.sub(
            r'^[ISHP]+', '',
            re.sub(
                r'[ISHP]+$', '',
                ''.join([
                    (k * int(v)) for k, v in zip(
                        re.sub(r'[0-9]', '', cigar),
                        re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
                    )
                ])
            )
        )
    )
    not_match_md = [
        (('D' * len(s[1:])) if s.startswith('^') else ('X' if s else s))
        for s in re.split(r'[0-9]+', md)
    ]
    n_not_match_md = len(not_match_md)
    match_md = [
        ('=' * int(s)) for s in re.split(r'[\^A-Z]', md) if s
    ]
    if n_not_match_md > len(match_md):
        match_md.append('')
    assert n_not_match_md == len(match_md)
    md_cigars = ''.join([(u + m) for m, u in zip(match_md, not_match_md)])
    ri = 0
    for i, c in enumerate(aln_cigar_list):
        if c == 'M':
            aln_cigar_list[i] = md_cigars[ri]
            ri += 1
        elif c in 'DN=X':
            ri += 1
    return ''.join(aln_cigar_list)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
