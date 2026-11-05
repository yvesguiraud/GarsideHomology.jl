# Resolutions, cells and elements

```@meta
CurrentModule = GarsideHomology
```

This page works in the reduced lcm resolution of the classical braid monoid B⁺(A₂), which has one cell of degree 0, two of degree 1 and one of degree 2. Its outputs are those of the tests of the package.

```julia
using GarsideHomology, Chevie

A2 = germ(coxgroup(:A, 2))       # simples e, 1, 2, 12, 21, 121
P = lcm_resolution(A2)
```

## Elements of a germ

A germ is called to give its elements, as [`GermElement`](@ref)s: `A2(1, 2)` is the product of the atoms 1 and 2, `A2("12")` the same simple read from its written form, `A2()` the identity, and `simples(A2)` the list of the simples. The product `s * t` is defined when it is a simple (`A2(1) * A2(2) == A2("12")`; `A2(1) * A2(1)` is an `ArgumentError`), and `s \ t` is the left quotient when `s ≼ t`. Also `word(s)`, the smallest word of `s` in the atoms, and `length(s)`.

## Cells

The cells of degree `n` are the basis of the module `P[n]`: `basis(P[1])` holds the cells `1` and `2`, and `basis(P)` all the cells, degree by degree. A cell is also given by its bar form, `x = P("2|12")`, or by its letters, `P(A2(2), A2(1, 2))`; `P()` is the cell of degree 0. A [`Cell`](@ref) reads as the sequence of its letters (`x[1]`, `x[end]`, `length(x)`), with its weight `cell_weight(x)`, here `121`, its degree `cell_degree(x)` and its name in the snapshots, `atom_sequence(x)`.

At the REPL, the resolution shows its cells and boundaries (abridged here):

```julia-repl
julia> P
LCMResolution: 4 cells, by degree [1, 2, 1], …
Cells:
  0: e
  1: 1, 2
  2: 2|12
Differential:
  (boundaries of degree 1 omitted, such as d(1) = 1[] - []1)
  d(2|12) = …
```

`show(stdout, "text/plain", P)`, or the same into a file, writes all the boundaries.

## Elements of the modules

The elements of the modules are sums of terms `c×σ[x]τ`: a cell `x`, an integer `c` and simples `σ` and `τ` acting on the two sides ([`BasedFreeModuleElement`](@ref)). They are built from cells, read from their written form, and combined:

```julia
y = P(P("1"))                    # the generator [1]
P("2×1[2]1 - [1]")               # written -[1] + 2×1[2]1
A2(1) * P(P("2")) == P("1[2]")   # the action on the left
y * A2(2) == P("[1]2")           # and on the right
```

The identity is omitted on both sides, the coefficients ±1 too, and the cell of degree 0 is written `[]`. The marker `×` separates a coefficient from what follows when the atoms are written with digits. `v[σ, x, τ]` is the coefficient of a term, and the iteration over `v` gives the terms as tuples `(c, σ, x, τ)`.

## The differential

```julia
d = differential(P)
d(P("1"))                        # 1[] - []1
d(A2(2) * y) == P("21[] - 2[]1")
iszero(d(d(x)))                  # true
```

`d(x)` is the boundary of the generator `x`, and `d(v)` that of an element ([`Differential`](@ref)). The matrices of the differential are `differential(P, n)`, from degree `n + 1` to degree `n`: [`LabelledSparseMatrix`](@ref)es whose entries carry the pairs `(σ, τ)`.

## Extensions of scalars

`E = extend_scalars(P)` is the complex over ℤ obtained by letting the simples act by 1 on both sides, `ℤ ⊗ P ⊗ ℤ`; its bases are the same cells, and its matrices are `SparseMatrixCSC`s. Its elements are written without simples (`E(P("1"))` is `[1]`), and `Vector(differential(E)(x))` gives the column of `x`. The other tensor products, with ℤε for instance, are on the page [Coefficients and homology](homology.md).
