# Getting started

```@meta
CurrentModule = GarsideHomology
```

This page goes from a reflection group to the homology of its braid group. The monoids are built with [Chevie.jl](https://github.com/jmichel7/Chevie.jl); the installation is described on the [Home](../index.md) page.

## A first computation

```julia
using GarsideHomology, Chevie

W = crg(4)                    # the complex reflection group G₄
G = dual_germ(W)              # the germ of its dual braid monoid B*(G₄)
R = lcm_resolution(G)         # its reduced lcm resolution
certify(R)                    # correct, minimal, linear
homology(R, ℤ)                # ℤ, ℤ, 0 in degrees 0 to 2
homology(R, ℤε)               # ℤ/2ℤ, ℤ/3ℤ, 0
```

The steps are the following.
- A *germ* gives the simples of a Garside monoid and their products. [`dual_germ`](@ref) builds the germ of the dual braid monoid of a well-generated irreducible reflection group, without building the monoid itself; [`germ`](@ref) builds that of a classical braid monoid, `germ(coxgroup(:A, 3))`, or of a monoid of Chevie, `germ(BraidMonoid(W))`.
- [`lcm_resolution`](@ref) builds a free resolution of the monoid, a chain complex of free bimodules whose bases are cells.
- [`certify`](@ref) checks the resolution (homogeneity of the boundaries, ∂² = 0, Euler characteristic of each weight) and says whether it is minimal.
- `homology(R, M)` is the homology of the braid group with coefficients in `M`: `ℤ` for the trivial coefficients, `ℤε` for the sign module, in which the braid reflections act by −1.

The result is a [`FinitelyGeneratedGradedAbelianGroup`](@ref), shown one degree per line:

```julia-repl
julia> homology(R, ℤ)
FinitelyGeneratedGradedAbelianGroup in degrees 0:2
  0: ℤ
  1: ℤ
  2: 0
```

`homology(R, ℤ, 1)` gives the degree 1 alone. The page [Background and conventions](background.md) recalls the notions and the conventions used here.

## Shortcuts

`homology(G, ℤε)` builds the resolution of the germ `G` and computes its homology in one call, and with Chevie loaded, `homology(BraidMonoid(W))` does the same for a monoid of Chevie. Each of these calls builds the resolution again: to compute several homologies of the same monoid, build the resolution once with `lcm_resolution`.

For the dual braid monoids, `dual_germ(W)` is much faster than `germ(DualBraidMonoid(W))`: building `DualBraidMonoid(W)` itself with Chevie takes minutes for the largest groups (about 4 for G₃₄), where `dual_germ(W)` takes a fraction of a second.

## Classical and dual braid monoids

- `germ(coxgroup(:E, 6))` is the germ of the classical braid monoid B⁺(E₆), whose simples are tabulated. For a Coxeter group with more than 5·10⁵ elements, such as E₈, `germ` gives a lazy [`CoxeterGerm`](@ref) instead, which numbers the simples as the construction meets them.
- `dual_germ(W)` accepts Coxeter groups (`coxgroup(:B, 4)`) and well-generated complex reflection groups (`crg(24)`, `crg(4, 4, 3)`). The Coxeter element is by default the product of the generators in their order, when it is regular; otherwise, `dual_germ` takes the first order of the generators that gives a regular element, and says which (for G₂₉: `[1, 2, 4, 3]`). Another Coxeter element is given as a word in the generators: `dual_germ(coxgroup(:A, 3); c = [2, 1, 3])`. There is no lazy dual germ: above 5·10⁵ simples, `dual_germ` throws an `ArgumentError`.

## Threads and times

Start Julia with several threads, `julia -t auto`: the check of ∂² = 0 in [`certify`](@ref) and the computation of the elementary divisors use them. The first call of a session includes the compilation of the code. The time and the memory grow quickly with the size of the monoid, mostly through the number of cells and the size of the differentials.
