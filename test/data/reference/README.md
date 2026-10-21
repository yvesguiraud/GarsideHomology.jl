# Reference fingerprints (J0)

Fingerprints of the reference snapshots computed in J0 with the preliminary
implementation (`gauss.jl` of the author's experiment repository), before the
construction was rewritten. The new construction must reproduce them exactly: strict
comparison on the atom sequences, bar forms as a second level of diagnosis.

- 344 files `<case>_<variant>.sha.csv`: 86 cases (34 classical
  monoids B⁺(W), 26 dual monoids B*(W) of real groups, 26 dual monoids of well-generated
  complex reflection groups), each in the four variants a to d.
- Format of the snapshots and of these files: [FORMAT.md](FORMAT.md) (version 1).
- The snapshots themselves (`.txt.gz`) are not in this repository: 72 MB for the 300 small
  ones, and more for the heavy ones. The comparison only needs these fingerprints.

Provenance, read from the headers of the 336 snapshots kept in the experiment
repository (the snapshots of G37, dual_E8, heavier, are kept outside git):

- Chevie: 0.1.15
- Julia: 1.13.1
- SHA-256 of `gauss.jl` (LF line endings): a45f0c3b02d217da3a0287de6fc767e57bf134b4f5f65a40d1ff2ce264fd0b9c
- Commit of the experiment repository: 6ea7f2ecc8cb36c40936b43a0e963e9871de9998, 8a2ba617063dc85663aca7990b57637c260a8f8d, cad379d4c3a7fe9898cdded1492208409ef2cb21
- Computed from 2026-10-14 to 2026-10-15.

## Collapsible pairs

`collapsible_pairs.csv` gives, for each snapshot, the counts of collapsible pairs recorded
in J0 by the experiment repository (script B3), extracted by its script
`j1/b4_collapsible_reference.jl`: the number of cells x such that ∂x has a term 1·y·1
(`sources`), the number of such pairs (x, y) (`pairs`) and the number of distinct cells y
(`targets`). Columns: `case,variant,sources,pairs,targets`.

The column `targets` is empty for the 248 snapshots recorded only by the first series of
J0, whose output did not count the targets. Where both series recorded a snapshot, they
agree on sources and pairs.

These terms were counted whatever their coefficient. The tests compare these counts with
the monomials of label (e, e) of the snapshots, and check that `collapsible_pairs` keeps
exactly those whose coefficient is ±1.
