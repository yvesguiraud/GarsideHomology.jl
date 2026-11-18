# Getting started

```@meta
CurrentModule = GarsideHomology
```

This page goes from a reflection group to the homology of its braid group. The monoids are built with [Chevie.jl](https://github.com/jmichel7/Chevie.jl); the installation is described on the [Home](../index.md) page.

## A first computation

Launch Julia as usual, preferably with several threads (`julia -t auto`) to speed up the computation of the elementary divisors (done by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl)) and the acyclicity check of complexes (d² = 0). As usual with Julia, the fist call of each function in a session includes the compilation of the code.

```julia
using GarsideHomology, Chevie

W = crg(4)                     # the complex reflection group G₄ (constructed by Chevie)
G = dual_germ(W)               # the germ of the dual braid monoid B*(G₄)
P = lcm_resolution(G)          # its reduced lcm resolution
certify(P)                     # correct, minimal, linear
homology(P, ℤ)                 # ℤ, ℤ, 0 in degrees 0 to 2
homology(P, ℤε)                # ℤ/2ℤ, ℤ/3ℤ, 0 in degrees 0 to 2
```

The steps are the following:
- the construction of a *germ* of the Garside monoid: a structure containing the simples, their products and other useful operations on them;
- the construction of the *lcm resolution* of the germ;
- the computation of the homology of the monoid in a given module, using the lcm resolution.

## Construction of the germ

The germ can be built using one of the two functions:
- [`germ`](@ref) builds the germ of a classical or dual braid monoid built by Chevie;
- [`dual_germ`](@ref) is an alternative construction, building the germ of the dual braid monoid of a finite Coxeter group or of a well-generated irreducible reflection group, without building the monoid itself.

For example, `germ(BraidMonoid(W))` works for every finite Coxeter group `W` built using the Chevie functions `coxgroup`, and `germ(DualBraidMonoid(W))` works for every finite Coxeter group and for every well-generated complex reflection group `W` built using the Chevie function `crg`. There is a shortcut `germ(W)` for the germ of the classical braid monoid B⁺(`W`): for example, `germ(coxgroup(:E, 6))` returns the germ of the classical braid monoid B⁺(E₆). For Coxeter monoids with few elements, `germ(W)` computes the entire germ of B⁺(`W`); over a certain threshold, empirically fixed at 5·10⁵ elements at the moment (thus including E₈ for example), `germ(W)` returns a lazy [`CoxeterGerm`](@ref) instead, which stores the simples and operations on them only when they are encountered in subsequent computations.

The alternative construction `dual_germ(W)` is much faster than `germ(DualBraidMonoid(W))`: building `DualBraidMonoid(W)` itself with Chevie takes minutes for the largest groups (about 4 minutes for G₃₄), whereas `dual_germ(W)` takes a fraction of a second. The function `dual_germ` accepts Coxeter groups (`W = coxgroup(:B, 4)`) and well-generated complex reflection groups (`W = crg(24)`, `W = crg(4, 4, 3)`). The Coxeter element is by default the product of the generators of `W` in the order fixed by Chevie, when this product is regular; otherwise, `dual_germ` takes the first order of the generators that gives a regular element, and displays which one (for example, it takes `[1, 2, 4, 3]` for G₂₉). A custom Coxeter element can be specified as a word in the generators: `dual_germ(coxgroup(:A, 3); c = [2, 1, 3])`. There is no lazy dual germ: above 5·10⁵ simples, `dual_germ` throws an `ArgumentError`.

## Construction of the lcm resolution

The main function [`lcm_resolution`](@ref) takes the germ of a Garside monoid M as argument and returns the reduced lcm resolution of M, which is a resolution of the monoid algebra ℤM by free ℤM-bimodules. An optional argument allows to specify an order on atoms and their lcms: the default value gives the reduced lcm resolution, but another order can be specified to recover the original order resolution of [DL99]. The package provides the shortcut `lcm_resolution(M)` where `M` can be a `BraidMonoid` or a `DualBraidMonoid`, which first computes `germ(M)`.

The test [`certify`](@ref) checks the resolution: homogeneity of the boundaries, d² = 0, Euler characteristic of each weight. It also says whether it is minimal, in the sense that no generating cells can be removed from the resolution.

## Computation of the homology

If `P` is the lcm resolution of a classical or dual braid monoid, then `homology(P, X)` gives the homology of the corresponding braid group with coefficients in the module `X`; for the moment, `X` can be either the trivial module `ℤ`, or the signature module `ℤε`, in which the braid reflections act by −1.

The result of `homology(P, X)` is a [`FinitelyGeneratedGradedAbelianGroup`](@ref), shown one degree per line:

```julia-repl
julia> homology(P, ℤ)
FinitelyGeneratedGradedAbelianGroup in degrees 0:2
  0: ℤ
  1: ℤ
  2: 0
```

If `H = homology(P, ℤ)`, then `H[n]` gives the `n`th homology group, stored as a [`FinitelyGeneratedAbelianGroup`](@ref), and `homology(P, ℤ, n)` computes it alone without computing the entire homology.

If `G` is a germ, then `homology(G, X)` is a shortcut that builds the resolution of the corresponding monoid and then computes its homology; `homology(M, X)` does the same for a classical or dual braid monoid built by Chevie. Each of these calls computes the whole resolution: to compute the homology of a given monoid in various coefficients, one should build the resolution once with `lcm_resolution` and apply `homology` on the result for each choice of coefficients.

The computation of `homology(P, X)` can be decomposed by first computing the complex `C = ℤ ⊗ P ⊗ X`, and then calling `homology(C)`.

## Typing the unicode characters

In the Julia REPL and in the editors that know Julia (VS Code, Jupyter, Pluto, the Julia modes of Emacs and Vim), unicode characters can be entered easily:
- `ℤ`: `\bbZ`, then Tab.
- `ℤε`: `\bbZ`, Tab, `\varepsilon`, Tab. For convenience, typing `ℤ` and `\epsilon` gives `ℤϵ`, which is recognised as the same module as `ℤε` by the package.
- `⊗`: `\otimes`, then Tab.