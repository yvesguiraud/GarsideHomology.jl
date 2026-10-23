# Format of the Gauss reference snapshots — version 1

Validated on 14 October 2026. Translated from the original French specification.

A snapshot fixes, for each case and each variant, the cells and the differential computed by
the preliminary implementation. The new construction must reproduce them. The format depends
only on the monoid, the variant and Chevie's order of the atoms: neither on the machine, nor on
the internal order of dictionaries, nor on line endings.

## 1. Files

| File | Contents |
| --- | --- |
| `<case>_<variant>.txt.gz` | the full snapshot, gzip-compressed (not in this repository) |
| `<case>_<variant>.sha.csv` | the fingerprints per degree, and that of the `.txt.gz` file |

`<case>` is an identifier such as `classic_E8`, `dual_B3` or `G34`, and `<variant>` a letter
from `a` to `d` (Gauss_a to Gauss_d).

## 2. Notation

- **Atoms**: indices 1, …, N in the order of `B.atoms` given by Chevie.
- **Simple**: its smallest reduced word in the atoms for the lexicographic order, indices
  separated by dots (`2.1.3`); the identity is written `e`.
- **Cell** x = [a₁|…|aₙ], with a₁⋯aⱼ = lcm(i₁, …, iⱼ), written in two ways:
  - *atom sequence* (i₁, …, iₙ): `7.3.12`; the cell of degree 0 is written `-`;
  - *bar form*: `[1.2|3|2.1]`; the cell of degree 0 is written `[]`.
- **Weight** π(x) = a₁⋯aₙ, written as a simple.
- **Term** of ∂x: coef · σ · y · τ, with σ the simple on the left and τ the simple on the
  right; dₙ : Cₙ₊₁ → Cₙ, so a term of ∂x belongs to dₙ when x has degree n + 1.
- **Coefficient**: signed integer, always with its sign (`+1`, `-1`, `+2`).

All data lines are ASCII; only the header contains UTF-8. LF line endings.

## 3. Header

Lines starting with `# `, of the form `# key: value`, in this order (keys are kept in French):

```
# format: cliche-gauss 1
# cas: dual_B3
# libelle: B*(B₃)
# type: dual                          (classic, dual or complex)
# construction: dual_braid_monoid(:B, 3)
# variante: a (filter = true, lex = true)
# cote: ppcm à droite ; terme coef·σ·y·τ, σ à gauche, τ à droite ; d_n : C_{n+1} -> C_n
# garside_W: 1.2.3                    (Coxeter element, or w₀ for a classical monoid, in generators of W)
# n_atomes: 9
# atome 1: 1                          (one line per atom, in generators of W)
# …
# atome 9: 2.1.2.3.2.1.2
# chevie: 0.1.15
# julia: 1.13.1
# commit_catlab: <commit of the experiment repository>
# sha256_gauss_jl_LF: <SHA-256 of gauss.jl with LF line endings>
# machine: <machine>
# date: 2026-10-14T16:21:23
# cellules: 38
# termes: 162
```

The comments in parentheses are explanations, not part of the files. The header does not
enter the fingerprints, since it contains the machine and the date.

## 4. Data lines

Cells, one line each:

```
C <degree> <atom sequence> <bar form> <weight>
```

Terms of the differential, one line each:

```
D <n> <source> <coef> <σ> <target> <τ>
```

where `<source>` and `<target>` are atom sequences. Schematic example (illustrative values):

```
C 0 - [] e
C 1 4 [4] 4
C 2 5.1 [5|1.3] 1.2.5
D 0 4 +1 4 - e
D 0 4 -1 e - 4
D 1 5.1 +1 e 5 1.3
```

Order of the file: the header, then the `C` lines by increasing degree, then the `D` lines by
increasing n. Within each group (same letter, same degree or same n), lines are sorted by bytes.

## 5. Fingerprints

For each group, the fingerprint is the SHA-256 of the concatenation of its sorted lines, each
followed by `\n`. It is computed in two forms:

- **atoms** (`sha256_atomes`): the lines as written in the file;
- **bar** (`sha256_barre`): the same lines where each cell is replaced by its bar form (for
  `C` lines, the atom-sequence field is removed), then sorted again.

The file `<case>_<variant>.sha.csv` has the columns

```
section,n,lignes,sha256_atomes,sha256_barre
```

with one line per group (`C` or `D`, degree or n, number of lines), one line `total` (numbers
of lines, and SHA-256 of the concatenation of the group fingerprints in file order), and one
line `fichier` giving the SHA-256 of the `.txt.gz` file.

## 6. Comparison

- **Strict**: all `sha256_atomes` fingerprints are equal.
- **Diagnosis**: if an `atoms` fingerprint differs but the `bar` fingerprints are equal, only
  the choice of the atom representing an lcm has changed. This is the case, for instance,
  between the tie-breaking rules of Gauss_a (first in the inherited order) and Gauss_c
  (smallest index).

## 7. Header of the snapshots written by Gauss.jl

Since J1.d (task d4, 23 October 2026), Gauss.jl writes snapshots itself
(`Gauss.write_snapshot`, for a resolution given by `lcm_resolution`). The data lines and
the fingerprints follow sections 4 and 5; the header, which does not enter the
fingerprints, differs from that of J0, which came from Chevie and from the experiment
repository. Keys, in this order:

```
# format: cliche-gauss 1
# cas: dual_B3                          (optional, given by the caller)
# ordre: LCMOrdering(ComplementDeglex())
# cote: ppcm à droite ; terme coef·σ·y·τ, σ à gauche, τ à droite ; d_n : C_{n+1} -> C_n
# n_atomes: 9
# empreinte_germe: GermFingerprint(9 atoms, ℓ(Δ) = 3, …)   ("non calculée (germe paresseux)" for a lazy germ)
# …                                     (the pairs given by the caller: libelle, construction, machine, commits, …)
# gauss: 0.1.0-DEV
# julia: 1.13.2
# date: 2026-10-23T15:00:00
# cellules: 38
# termes: 162
```

The orderings that correspond to the variants of J0: Gauss_a is
`LCMOrdering(ComplementDeglex())` (identical to J0 except for the classical monoids from
E�� on, whose children are ordered by deglex instead of lex), Gauss_d is
`LCMOrdering(TrivialPreorder())`; Gauss_b and Gauss_c exist only as orderings of the tests.
