---
name: align-sequences
description: Align protein sequences and inspect the result as a table with pdbframe's SeqFrame - pairwise or several sequences, identity/similarity, gaps, finding a sequence inside structure chains. Use for sequence comparison questions.
---

# Align sequences with SeqFrame

```python
from pdbframe import SeqFrame
out = SeqFrame.align(seq1, seq2)                 # or SeqFrame(seq1, name='hla').align(seq2)
out = SeqFrame.align(s1, s2, s3, names=['a', 'b', 'c'])
out.df            # one row per column: residues, 1-based positions (<NA> at gaps), match
out.statistics    # identity, similarity, coverage, mismatches, gap_columns, score, per-seq start/end
out.show()        # text blocks; out.show(one_line=True)
out.drop_gaps()   # columns where every sequence has a residue; drop_gaps('a') for one
out.show(reference='a')   # others print '.' where identical to a
out.pairwise()    # {name: SeqFrame}: each sequence vs the first (or pairwise('a')), own statistics
```
`match` values: identical / similar / mismatch / gap. `mode='local'` (two sequences only)
keeps just the aligned region - use it to find a fragment inside a longer chain.

**From structures**: `pf.sequences` is `{chain: sequence}`; `pf.seqframe('A')` gives a SeqFrame,
`pf.seqframes()` one per chain. FASTA: `PDBFrame.read('x.pdb').write('x.fasta')`.

**Limits**
- Three or more sequences use a star alignment on the first one: fine for variants of one
  protein, not for distant families - use MAFFT there (`module load mafft` on the cluster).
- Identity is over alignment columns, so end gaps lower it; check `coverage` too.
- `drop_gaps()` is a column filter, not a local alignment: it also removes internal indels,
  and its `coverage` is over the remaining residues. For identity against the first
  sequence per pair in a multiple alignment, use `pairwise()`.
- Scoring defaults: match 2, mismatch -1, gap open -2, gap extend -1 (keyword arguments).

Every method and option: `pdbframe docs seqframe` (prints `docs/seqframe.md`).
