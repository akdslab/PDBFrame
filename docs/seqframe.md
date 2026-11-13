# SeqFrame reference

A `SeqFrame` is a sequence, or an alignment of several sequences, as a pandas table with
one row per position. Every method returns a new SeqFrame (or a table); nothing is changed
in place.

```python
from pdbframe import PDBFrame, SeqFrame

pf = PDBFrame.from_rcsb('3tfk')                          # an mmCIF file
aln = pf.seqframe('A', full=True).align(pf.seqframe('A'), names=['deposited', 'modeled'])
aln.show()
```

`full=True` takes the full sequence, unobserved residues included (`pf.full_sequences`);
without it, the residues with atoms (`pf.sequences`). A PDB file has no full sequence; a
FASTA file (`PDBFrame.read('3TFK.fasta')`) gives the deposited sequences on their own.

```
deposited     1 MGPHSMRYYETATSRRGLGEPRYTSVGYVDDKEFVRFDSDAENPRYEPQVPWMEQEGPEY 60
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
modeled       1 --PHSMRYYETATSRRGLGEPRYTSVGYVDDKEFVRFDSDAENPRYEPQVPWMEQEGPEY 58
...
deposited   121 GCDYIALNEDLRTWTAADMAAQITRRKWEQAGAAEYYRAYLEGECVEWLHRYLKNGNATL 180
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
modeled     119 GCDYIALNEDLRTWTAADMAAQITRRKWEQAGAAEYYRAYLEGECVEWLHRYLKNG---- 174
```

## Building one

| | |
|---|---|
| `SeqFrame('MKVL', name=None)` | one sequence; `name` is its column name (default `'seq1'`) |
| `pf.seqframe('A', full=False)` | one chain of a PDBFrame, named after the chain (`'A'`); `full=True` with its unobserved residues |
| `pf.seqframes(full=False)` | every chain of a PDBFrame, `{chain: SeqFrame}` (or `pf.seqframes(['A', 'B'])`) |
| `SeqFrame.ensure(x)` | a string or a SeqFrame, as a SeqFrame |

## The table

`sf.df` for one sequence has a single column of residues, indexed by position `pos`
(1, 2, 3, ...).

An alignment is indexed by alignment column `aln_pos` and has, for every sequence:

| column | |
|---|---|
| `<name>` | the residue, or `-` across a gap |
| `<name>_pos` | 1-based position of that residue in its own sequence; `<NA>` across a gap |
| `match` | `'identical'`, `'similar'` (same residue group, e.g. D/E or K/R/H), `'mismatch'` or `'gap'` |

## Aligning

| | gives |
|---|---|
| `SeqFrame.align(seq1, seq2, ...)` or `sf.align(other, ...)` | the alignment; the two calls are the same |

Two sequences get a Biopython pairwise alignment. Three or more get a star alignment: each
is aligned to the first and the results are merged, which is fast and suits related
sequences (it is not a full multiple aligner like MAFFT or Clustal).

Options:
- `mode='global'` (default) or `'local'`, which keeps only the aligned region (two
  sequences only).
- `names=[...]`: column names for all sequences, in order. Defaults to each SeqFrame's
  name, or `seq1`, `seq2`, ...
- `match_score=2`, `mismatch_score=-1`, `gap_open=-2`, `gap_extend=-1`: the scoring.

```python
SeqFrame.align('CASSIRSSYEQYF', 'CASSIRSAYEQYF', 'CASSLGQAYEQYF').show()
```

```
seq1     1 CASSIRSSYEQYF 13
seq2     1 CASSIRSAYEQYF 13
seq3     1 CASSLGQAYEQYF 13
           ||||:...|||||
```

### Against one reference

With three or more sequences the first is the reference, and every row above is placed
relative to it. Two more ways to look at that:

```python
cdr3 = SeqFrame.align('CASSIRSSYEQYF', 'CASSIRSAYEQYF', 'CASSLGQAYEQYF', 'CASSIRSSGYEQYF',
                      names=['tcr1', 'tcr2', 'tcr3', 'tcr4'])
cdr3.show(reference='tcr1')      # residues identical to tcr1 print as '.'
```

```
tcr1     1 CASSIRSS-YEQYF 13
tcr2     1 .......A-..... 13
tcr3     1 ....LGQA-..... 13
tcr4     1 ........G..... 14
           ||||:... |||||
```

`cdr3.pairwise()` gives each sequence against the reference as its own two-sequence
alignment, `{name: SeqFrame}`, so each has its own match line and statistics:

```python
for p in cdr3.pairwise().values():
    p.show()
pd.DataFrame({n: p.statistics for n, p in cdr3.pairwise().items()})   # one column per pair
```

```
tcr1     1 CASSIRSSYEQYF 13
           |||||||.|||||
tcr2     1 CASSIRSAYEQYF 13
...
tcr1     1 CASSIRSS-YEQYF 13
           |||||||| |||||
tcr4     1 CASSIRSSGYEQYF 14
```

Each pair keeps the columns where either sequence has a residue, so the gap that tcr4's
`G` opens in the multiple alignment is gone from the tcr1-tcr2 pair. For an alignment
from `align`, the pairs against the first sequence are the pairwise alignments it was
built from. `pairwise('tcr3')` uses another sequence as the reference (its pairs are cut
from this alignment, not re-aligned).

## Contents

| | gives |
|---|---|
| `sf.names` | the sequence column names |
| `sf.sequences` | `{name: sequence}` without gaps |
| `sf.aligned` | `{name: sequence}` with gaps |
| `sf.is_alignment` | True when it holds more than one sequence |
| `sf.score` | the alignment score (pairwise only, else `None`) |
| `sf.pairwise(reference=None)` | `{name: SeqFrame}`, each sequence against the reference (the first by default) |
| `len(sf)` | number of positions or alignment columns |

## Statistics

`sf.statistics` (a property, not a method) is a Series:
- one sequence: name, type (`'protein'`, `'dna'` or `'rna'`), length and residue counts.
- an alignment: `identity` and `similarity` (share of columns that are all identical /
  all similar), `coverage` (identical columns over the first sequence's length), the
  counts of matches, similar, mismatches and gap columns, the score, and each sequence's
  gap count and start / end.

For the 3TFK example above, identity is 0.97: 174 of the 180 deposited residues are in
the structure.

## Editing

| | gives |
|---|---|
| `sf.drop_gaps()` | only the columns where every sequence has a residue |
| `sf.drop_gaps('modeled')` | only the columns where `modeled` has a residue |
| `sf.drop_gaps(how='all')` | drops a column only if every sequence has a gap there |

The `aln_pos` index is kept, so you can see which columns survived. The score is cleared.

`drop_gaps()` is not a local alignment. It filters the columns of the alignment you
already have; it does not re-align, so `statistics['mode']` stays `'global'`. It drops
*every* gapped column, the gaps inside (insertions and deletions) as well as the
overhangs at the ends. `mode='local'` re-aligns to find the best-matching region, and
that region can keep gaps inside it. The two give the same table only when all the gaps
are at the ends, as in the notebook example.

So the statistics of `drop_gaps()` describe the gap-free columns only: `identity` is the
identity over the residues the sequences share. `coverage` is computed on the table as
it stands, so it is over the first sequence's remaining residues, not its full length;
take coverage from the alignment before dropping gaps.

## Printing

| | |
|---|---|
| `sf.show(width=60)` | print in blocks of `width`, with start / end positions. The match line shows `\|` identical, `:` similar, `.` mismatch |
| `sf.show(one_line=True)` | everything on one line |
| `sf.show(reference='tcr1')` | the other sequences print `.` where they are identical to `tcr1` (`True`: the first sequence) |
| `sf` in a notebook | the table |
