# Homology

```@meta
CurrentModule = GarsideHomology
```

## Modules of coefficients

A module of coefficients is a module over the monoid algebra, free of finite rank over ``ℤ`` ([`AbstractMonoidModule`](@ref)). Two are defined for now:
- [`TrivialModule`](@ref)`()`, or `ℤ` for short, in which every element of the monoid acts trivially;
- [`SignModule`](@ref)`()`, or `ℤε` or `ℤϵ` for short, in which every atom acts by ``−1`` (for a braid monoid, every braid reflection), so that a simple `s` acts by `(-1)^length(s)`.

## Homology

If `P` is a resolution of a Garside monoid ``M`` and `X` is a module of coefficients, then `homology(P, X)` computes the homology of ``M`` with coefficients in `X`, in all degrees, and `homology(P, X, n)` computes only the homology in degree `n`. Since ``M`` is a Garside monoid, if `X` can be extended to a module over the enveloping group ``G(M)`` of ``M`` (which is the case for the trivial and sign modules), the computed homology is also that of ``G(M)`` with coefficients in `X` (see [Mathematical background](background.md)).

The result `H` of `homology(P, X)` is a [`FinitelyGeneratedGradedAbelianGroup`](@ref). The individual homology groups can be accessed through `H[n]`, which is a [`FinitelyGeneratedAbelianGroup`](@ref). The free rank and the torsion of `H[n]` are given by `rank(H[n])` and `torsion(H[n])`, respectively. Finitely generated abelian groups can be entered using the `parse` function, as in `parse(FinitelyGeneratedAbelianGroup, "(ℤ/3ℤ)² ⊕ ℤ/6ℤ ⊕ ℤ")`; the function accepts `Z` for `ℤ`, `ℤ/2` for `ℤ/2ℤ`, `^k` for an exponent and `+` for `⊕`.

The elementary divisors of the differentials are computed by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl), and the keyword arguments of `homology` are transferred to it, such as `value_type`, for the integer types tried during the elimination, `nworkers`, for the number of workers to be used (if multithreading is activated), `verbose`, and the others described in its documentation.

## Tensor products

The function `homology(P, X)` returns the homology of the complex `ℤ ⊗ P ⊗ X` of abelian groups. The function [`tensor_product`](@ref), or `⊗` for short, can be used to compute other complexes:

```julia
using GarsideHomology: ⊗           # not necessary if Chevie is loaded

C = ℤ ⊗ P ⊗ ℤε                     # a complex of abelian groups
L = P ⊗ ℤ                          # a resolution by free left modules
N = ℤ ⊗ P                          # a resolution by free right modules
```

The special case `C = ℤ ⊗ P ⊗ ℤ` can be built directly by `C = extend_scalars(P)`: it is the complex of abelian groups obtained by letting the simples act trivially on both sides. The basis of `C` is the same as that of `P`, and the matrices of its differential are `SparseMatrixCSC`s. The elements of the modules composing `C` are written without simples (e.g. `C("1")` is `[1]`), and `Vector(differential(C)(x))` gives the column corresponding to `x`.

## Complexes of abelian groups

The generic layer, `GarsideHomology.Complexes`, also computes the homology of complexes of integer matrices, given by their differentials `d_0, d_1, …` with `d_n : C_{n+1} → C_n`, and not only of those coming from the lcm resolution:

```julia
using GarsideHomology, SparseArrays

# Cellular chain complex of the real projective plane: d₀ = 0, d₁ = 2.
C = ChainComplex([spzeros(Int, 1, 1), sparse([1], [1], [2], 1, 1)])
homology(C)                        # ℤ, ℤ/2ℤ, 0 in degrees 0 to 2
is_square_zero(C)                  # true
```

The function [`collapsible_pairs`](@ref) lists the pairs of basis elements joined by a coefficient ``±1`` in a differential, for a complex of abelian groups or of (left, right, bi)-modules over a monoid algebra.