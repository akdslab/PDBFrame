"""pandas-based sequence object.

A :class:`SeqFrame` is a per-position table. A single sequence has one column of
residues. An alignment has one column per sequence, each sequence's 1-based residue
numbers, and a ``match`` label per column::

    >>> out = SeqFrame('MKVLWAAL').align('MKVAAL')      # or SeqFrame.align(s1, s2)
    >>> out                    # the alignment as a DataFrame
    >>> out.statistics         # identity, similarity, gaps, score, ...
    >>> out.show()             # printed alignment; out.show(one_line=True) for one line
    >>> out.drop_gaps()        # only the columns where every sequence has a residue
    >>> out.pairwise()         # each sequence against the first, as pairwise alignments

More than two sequences align the same way, ``SeqFrame.align(s1, s2, s3)``.
"""

from __future__ import annotations

import functools
from collections import Counter

import pandas as pd
from Bio import Align

__all__ = ['SeqFrame', 'detect_type']

GAP = '-'

# Residues that count as "similar" (':') when proteins disagree.
SIMILARITY_GROUPS = [set('GAVLI'), set('FYW'), set('CM'), set('ST'), set('NQ'),
                     set('DE'), set('KRH')]

# Values of the match column, and the symbol show() draws for each.
MATCH, SIMILAR, MISMATCH, GAPPED = 'identical', 'similar', 'mismatch', 'gap'
SYMBOLS = {MATCH: '|', SIMILAR: ':', MISMATCH: '.', GAPPED: ' '}


def detect_type(sequence: str) -> str:
    """Guess ``'protein'``, ``'dna'`` or ``'rna'`` from the characters present."""
    chars = set(sequence)
    if chars & set('EFILPQZ'):
        return 'protein'
    if 'U' in chars:
        return 'rna'
    if 'T' in chars:
        return 'dna'
    return 'protein'


def _label(residues: tuple[str, ...], protein: bool) -> str:
    """The match label for one alignment column."""
    chars = set(residues)
    if GAP in chars:
        return GAPPED
    if len(chars) == 1:
        return MATCH
    if protein and any(chars <= g for g in SIMILARITY_GROUPS):
        return SIMILAR
    return MISMATCH


def _aligned_strings(alignment, s1: str, s2: str) -> tuple[str, str]:
    """The two gapped rows, rebuilt from the alignment coordinates.

    Done by hand rather than with ``str(alignment[0])`` because that call returns
    different spans for local alignments across Biopython versions.
    """
    row1, row2 = [], []
    coords = alignment.coordinates
    for (a0, b0), (a1, b1) in zip(coords.T[:-1], coords.T[1:]):
        row1.append(s1[a0:a1] if a1 > a0 else GAP * (b1 - b0))
        row2.append(s2[b0:b1] if b1 > b0 else GAP * (a1 - a0))
    return ''.join(row1), ''.join(row2)


def _merge_star(ref: str, pairs: list[tuple[str, str]]) -> list[str]:
    """Merge pairwise global alignments to `ref` into one multiple alignment.

    Each pair is (gapped ref, gapped other). Residues opposite a ref residue stay
    there; residues inserted between two ref residues are stacked into the widest
    insertion any sequence needs at that spot, padded with gaps.
    """
    L = len(ref)
    inserts, opposite = [], []          # per other: chars before ref[j], char at ref[j]
    for r, o in pairs:
        ins, opp, j = [[] for _ in range(L + 1)], [], 0
        for a, b in zip(r, o):
            if a == GAP:
                ins[j].append(b)
            else:
                opp.append(b)
                j += 1
        inserts.append(ins)
        opposite.append(opp)

    rows = [[] for _ in range(len(pairs) + 1)]
    for j in range(L + 1):
        width = max(len(ins[j]) for ins in inserts)
        if width:
            rows[0].append(GAP * width)
            for k, ins in enumerate(inserts, 1):
                rows[k].append(''.join(ins[j]).ljust(width, GAP))
        if j < L:
            rows[0].append(ref[j])
            for k, opp in enumerate(opposite, 1):
                rows[k].append(opp[j])
    return [''.join(r) for r in rows]


def _positions(aligned: str, start: int) -> pd.array:
    """1-based residue number at each alignment column, ``<NA>`` across a gap."""
    out, i = [], start
    for c in aligned:
        if c == GAP:
            out.append(pd.NA)
        else:
            i += 1
            out.append(i)
    return pd.array(out, dtype='Int64')


def _unique(names: list[str]) -> list[str]:
    """Suffix repeated names ``_1``, ``_2``... so every column is distinct."""
    counts = Counter(names)
    seen = Counter()
    out = []
    for n in names:
        seen[n] += 1
        out.append(f'{n}_{seen[n]}' if counts[n] > 1 else n)
    return out


class _class_or_instance:
    """Lets ``align`` be called as ``SeqFrame(s1).align(s2)`` or ``SeqFrame.align(s1, s2)``."""

    def __init__(self, func):
        self.func = func
        functools.update_wrapper(self, func)

    def __get__(self, obj, cls):
        if obj is None:
            @functools.wraps(self.func)
            def bound(first, *others, **kwargs):
                return self.func(cls.ensure(first), *others, **kwargs)
            return bound
        return functools.partial(self.func, obj)


class SeqFrame:
    """A sequence, or an alignment, as a per-position DataFrame.

    Args:
        sequence: a sequence string (or another SeqFrame, which is copied).
        name: column name for the sequence. Defaults to ``'seq1'``.

    Attributes:
        df: the table. Index is the 1-based position (``pos``) or alignment column
            (``aln_pos``).
        names: the sequence column names.
        score: the alignment score (pairwise only; ``None`` otherwise).
    """

    def __init__(self, sequence, name: str | None = None):
        if isinstance(sequence, SeqFrame):
            self.df, self.names = sequence.df.copy(), list(sequence.names)
            self.mode, self.score = sequence.mode, sequence.score
            return
        seq = str(sequence).upper().strip()
        name = name or 'seq1'
        self.df = pd.DataFrame({name: list(seq)},
                               index=pd.RangeIndex(1, len(seq) + 1, name='pos'))
        self.names = [name]
        self.mode = None
        self.score = None

    @classmethod
    def ensure(cls, data, name: str | None = None) -> 'SeqFrame':
        return data if isinstance(data, SeqFrame) else cls(data, name=name)

    @classmethod
    def _from_rows(cls, names, rows, starts, mode, score=None) -> 'SeqFrame':
        protein = all(detect_type(r.replace(GAP, '')) == 'protein' for r in rows)
        data = {n: list(r) for n, r in zip(names, rows)}
        data |= {f'{n}_pos': _positions(r, s) for n, r, s in zip(names, rows, starts)}
        data['match'] = [_label(col, protein) for col in zip(*rows)]
        out = cls.__new__(cls)
        out.df = pd.DataFrame(data, index=pd.RangeIndex(1, len(rows[0]) + 1, name='aln_pos'))
        out.names, out.mode, out.score = list(names), mode, score
        return out

    @property
    def is_alignment(self) -> bool:
        return len(self.names) > 1

    @property
    def sequences(self) -> dict[str, str]:
        """``{name: sequence}`` without gaps."""
        return {n: s.replace(GAP, '') for n, s in self.aligned.items()}

    @property
    def aligned(self) -> dict[str, str]:
        """``{name: sequence}`` as aligned, gaps included."""
        return {n: ''.join(self.df[n]) for n in self.names}

    def __len__(self) -> int:
        return len(self.df)

    def __str__(self) -> str:
        return self.df.to_string()

    def __repr__(self) -> str:
        return self.df.__repr__()

    def _repr_html_(self) -> str:
        return self.df._repr_html_()

    # -------------------------------------------------------------- alignment

    @_class_or_instance
    def align(self, *others, mode: str = 'global', names: list[str] | None = None,
              match_score: float = 2, mismatch_score: float = -1,
              gap_open: float = -2, gap_extend: float = -1) -> 'SeqFrame':
        """Align this sequence with one or more others.

        Two sequences get a Biopython pairwise alignment. Three or more get a star
        alignment: each is aligned to this (the first) sequence and the results are
        merged, so the first sequence acts as the reference. That is fast and
        dependency-free, and good for related sequences; it is not a full
        progressive aligner like MAFFT or Clustal.

        Args:
            others: sequence strings or SeqFrames.
            mode: ``'global'`` or ``'local'`` (only the aligned region is kept;
                two sequences only).
            names: column names for all sequences, in order. Defaults to the
                SeqFrames' own names, or ``seq1``, ``seq2``, ...
            match_score, mismatch_score, gap_open, gap_extend: scoring for
                Biopython's PairwiseAligner.

        Returns:
            A new SeqFrame whose table is the alignment: one row per alignment
            column, each sequence's residue and 1-based position (``<NA>`` across a
            gap), and the ``match`` label: ``'identical'`` (all residues the same),
            ``'similar'`` (same residue group), ``'mismatch'`` or ``'gap'``.
        """
        if not others:
            raise ValueError('align needs at least one other sequence')
        frames = [self] + [o if isinstance(o, SeqFrame) else SeqFrame(o, name=f'seq{i}')
                           for i, o in enumerate(others, 2)]
        if any(f.is_alignment for f in frames):
            raise ValueError('align takes single sequences, not an alignment')
        if mode == 'local' and len(frames) > 2:
            raise ValueError("mode='local' supports two sequences only")
        if names is not None and len(names) != len(frames):
            raise ValueError(f'{len(names)} names for {len(frames)} sequences')
        names = _unique(list(names) if names else [f.names[0] for f in frames])
        seqs = [next(iter(f.sequences.values())) for f in frames]

        aligner = Align.PairwiseAligner()
        aligner.mode = mode
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend
        # next() takes the top alignment without enumerating all co-optimal ones,
        # which can be astronomically many
        alignments = [next(iter(aligner.align(seqs[0], s))) for s in seqs[1:]]

        if len(frames) == 2:
            aln = alignments[0]
            rows = _aligned_strings(aln, seqs[0], seqs[1])
            starts = int(aln.coordinates[0][0]), int(aln.coordinates[1][0])
            return SeqFrame._from_rows(names, rows, starts, mode, float(aln.score))

        pairs = [_aligned_strings(a, seqs[0], s) for a, s in zip(alignments, seqs[1:])]
        rows = _merge_star(seqs[0], pairs)
        return SeqFrame._from_rows(names, rows, [0] * len(rows), mode)

    def drop_gaps(self, subset: list[str] | str | None = None,
                  how: str = 'any') -> 'SeqFrame':
        """Drop alignment columns that contain gaps.

        Args:
            subset: only look for gaps in these sequences, e.g. ``'seq1'`` keeps every
                column where seq1 has a residue. Defaults to all sequences.
            how: ``'any'`` drops a column if any of them has a gap, ``'all'`` only if
                all of them do.

        Returns:
            A new SeqFrame. The original ``aln_pos`` index is kept, so you can see
            which columns survived. The score is cleared, since it no longer
            describes the table.
        """
        cols = [subset] if isinstance(subset, str) else list(subset or self.names)
        gaps = self.df[cols] == GAP
        drop = gaps.any(axis=1) if how == 'any' else gaps.all(axis=1)
        out = SeqFrame(self)
        out.df = self.df[~drop]
        out.score = None
        return out

    def pairwise(self, reference: str | None = None) -> dict[str, 'SeqFrame']:
        """Each sequence against a reference, as two-sequence alignments.

        Every other sequence is paired with `reference` as they stand in this
        alignment; columns where both have a gap are dropped. For an alignment from
        :meth:`align` with the first sequence as the reference, these are the
        pairwise alignments it was built from, so each one's :attr:`statistics`
        give identity, coverage etc. against the reference alone.

        Args:
            reference: name of the reference sequence. Defaults to the first.

        Returns:
            ``{name: SeqFrame}``, one per other sequence, reference first in each.
            The score is ``None`` except for a two-sequence alignment, which is
            returned as a copy of itself.
        """
        if not self.is_alignment:
            raise ValueError('pairwise needs an alignment, not a single sequence')
        ref = self._reference(reference)
        out = {}
        for n in self.names:
            if n == ref:
                continue
            both_gap = (self.df[ref] == GAP) & (self.df[n] == GAP)
            cols = self.df[~both_gap]
            rows = [''.join(cols[ref]), ''.join(cols[n])]
            starts = [int(p.iloc[0]) - 1 if len(p := cols[f'{m}_pos'].dropna()) else 0
                      for m in (ref, n)]
            score = self.score if len(self.names) == 2 else None
            out[n] = type(self)._from_rows([ref, n], rows, starts, self.mode, score)
        return out

    def _reference(self, reference) -> str:
        if reference is None or reference is True:
            return self.names[0]
        if reference not in self.names:
            raise KeyError(f'no sequence {reference!r}; names are {self.names}')
        return reference

    # ------------------------------------------------------------- statistics

    @property
    def statistics(self) -> pd.Series:
        """Summary numbers.

        For one sequence: name, type, length and residue counts. For an alignment:
        identity / similarity (share of columns that are all identical / all
        similar), coverage (identical columns / length of the first sequence),
        column counts, the score (pairwise only) and each sequence's gap count and
        1-based start / end.
        """
        if not self.is_alignment:
            seq = self.sequences[self.names[0]]
            counts = dict(sorted(Counter(seq).items()))
            return pd.Series({'name': self.names[0], 'type': detect_type(seq),
                              'length': len(seq), **counts}, name='statistics')

        m = Counter(self.df['match'])
        width = len(self.df)
        first = self.sequences[self.names[0]]
        stats = {
            'sequences': len(self.names),
            'mode': self.mode,
            'score': self.score,
            'length': width,
            'identity': m[MATCH] / width if width else 0.0,
            'similarity': (m[MATCH] + m[SIMILAR]) / width if width else 0.0,
            'coverage': m[MATCH] / len(first) if first else 0.0,
            'matches': m[MATCH],
            'similar': m[SIMILAR],
            # every column with residues that are not all identical, similar included
            'mismatches': m[MISMATCH] + m[SIMILAR],
            'gap_columns': m[GAPPED],
        }
        for n in self.names:
            pos = self.df[f'{n}_pos'].dropna()
            stats[f'{n}_gaps'] = int((self.df[n] == GAP).sum())
            stats[f'{n}_start'] = int(pos.iloc[0]) if len(pos) else pd.NA
            stats[f'{n}_end'] = int(pos.iloc[-1]) if len(pos) else pd.NA
        return pd.Series(stats, name='statistics')

    # ----------------------------------------------------------------- display

    def show(self, width: int = 60, one_line: bool = False,
             reference: str | bool | None = None) -> None:
        """Print the sequence or alignment in blocks of `width`, or on one line.

        `reference` (a sequence name, or ``True`` for the first) prints the other
        sequences' residues that are identical to the reference's as ``.``, so only
        the differences show.
        """
        if one_line:
            width = max(len(self), 1)
        ref = None
        if reference is not None and reference is not False and self.is_alignment:
            ref = self._reference(reference)
        label = max(len(n) for n in self.names)
        # pairwise: match line between the two; more: match line underneath
        if len(self.names) == 2:
            rows = [self.names[0], 'match', self.names[1]]
        else:
            rows = self.names + (['match'] if self.is_alignment else [])
        for i in range(0, len(self), width):
            block = self.df.iloc[i:i + width]
            for r in rows:
                if r == 'match':
                    text = ''.join(SYMBOLS[v] for v in block[r])
                    print(f'{"":<{label}} {"":>5} {text}')
                    continue
                text = ''.join(block[r])
                if ref is not None and r != ref:
                    text = ''.join('.' if c == q and c != GAP else c
                                   for c, q in zip(block[r], block[ref]))
                pos = (block[f'{r}_pos'].dropna() if self.is_alignment
                       else block.index.to_series()).to_list()
                start, end = (pos[0], pos[-1]) if pos else ('', '')
                print(f'{r:<{label}} {start:>5} {text} {end}')
            print()
