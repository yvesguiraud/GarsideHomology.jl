# Tiny reference snapshots (J0)

Three of the reference snapshots of J0, small enough to be kept in the repository so
that the tests of the snapshot tool (`test/snapshot_io.jl`) run on real data in the CI:

| File | Monoid | Lines | Size |
| --- | --- | --- | --- |
| `classic_A2_a.txt.gz` | B⁺(A₂), Gauss_a | 14 | 519 bytes |
| `dual_A2_a.txt.gz` | B*(A₂), Gauss_a | 20 | 547 bytes |
| `dual_B3_a.txt.gz` | B*(B₃), Gauss_a | 200 | 1444 bytes |

They are copies of the files of the experiment repository from which the fingerprints
of `../reference/` were computed: their sizes and SHA-256 are those of the rows
`fichier` of the corresponding `.sha.csv`, which the tests check. Format:
[`../reference/FORMAT.md`](../reference/FORMAT.md).

The other snapshots are not in this repository. To check them, set the environment
variable `GAUSS_SNAPSHOT_DIR` to a folder of `.txt.gz` snapshots before running the tests.
