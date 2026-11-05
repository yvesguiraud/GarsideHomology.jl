# Coefficients and homology

```@meta
CurrentModule = GarsideHomology
```

## Modules of coefficients

A module of coefficients is a module over the monoid algebra, free of finite rank over ℤ ([`AbstractMonoidModule`](@ref)). Two are defined for now:
- [`TrivialModule`](@ref)`()`, in which every element of the monoid acts by 1;
- [`SignModule`](@ref)`()`, also written `ℤε`, in which every atom acts by −1, so that a simple `s` acts by `(−1)^ℓ(s)`.

A ring stands for its trivial module: `homology(R, ℤ)` is `homology(R, TrivialModule())`.

## Homology

For a resolution `R`, `homology(R, M)` is the homology of the monoid, hence of its group of fractions, with coefficients in `M`, in all degrees; `homology(R, M, n)` in degree `n`. The result is a [`FinitelyGeneratedGradedAbelianGroup`](@ref): `H[n]` is the [`FinitelyGeneratedAbelianGroup`](@ref) of degree `n`, with its free rank `rank(H[n])` and its invariant factors `torsion(H[n])`. A group is written torsion first, `(ℤ/3ℤ)² ⊕ ℤ/6ℤ ⊕ ℤ`, and this written form is read back by `parse(FinitelyGeneratedAbelianGroup, s)`, which also accepts `Z` for `ℤ`, `ℤ/2` for `ℤ/2ℤ`, `^k` for an exponent and `+` for `⊕`.

The elementary divisors are computed by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl), and the keyword arguments of `homology` go to it: `value_type`, the integer types tried during the elimination, `nworkers`, `verbose`, and the others described in its documentation.

## Tensor products

`homology(R, M)` is the homology of the complex `ℤ ⊗ R ⊗ M` over ℤ, which [`tensor_product`](@ref) and `⊗` give directly:

```julia
using GarsideHomology: ⊗           # or with Chevie loaded

C = ℤ ⊗ R ⊗ ℤε                     # the complex over ℤ
L = R ⊗ ℤ                          # free left modules: the terms keep their σ
N = ℤ ⊗ R                          # free right modules: the terms keep their τ
```

`ℤ ⊗ R ⊗ ℤ` is [`extend_scalars`](@ref)`(R)`. The results share the bases of `R`, whose elements are still cells, so that their elements are written and read as in `R`.

## Complexes over ℤ

The generic layer, `GarsideHomology.Complexes`, also computes the homology of complexes of integer matrices, given by their differentials `d_0, d_1, …` with `d_n : C_{n+1} → C_n`:

```julia
using GarsideHomology, SparseArrays

# Cellular chain complex of the real projective plane: d₀ = 0, d₁ = 2.
C = ChainComplex([spzeros(Int, 1, 1), sparse([1], [1], [2], 1, 1)])
homology(C)                        # ℤ, ℤ/2ℤ, 0 in degrees 0 to 2
is_square_zero(C)                  # true
```

[`collapsible_pairs`](@ref) lists the pairs of basis elements joined by a coefficient ±1 in a differential, for a complex over ℤ or over a monoid algebra.
