# Gauss.jl

[![CI](https://github.com/yvesguiraud/Gauss.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yvesguiraud/Gauss.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/yvesguiraud/Gauss.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yvesguiraud/Gauss.jl)

Small free resolutions of Garside monoids and their homology, in Julia.

Gauss.jl computes the *reduced lcm resolution* of a Garside monoid, a variant of the
Dehornoy–Lafont resolution with fewer cells, and the homology of the corresponding
groups with various coefficients. The first target is the homology of real and well-generated complex braid groups.

## Status

Early development, not registered. The first layer is in place: `Gauss.Complexes`, with
based free modules over monoid algebras, labelled sparse matrices, chain complexes, their
specialisation to integer coefficients and their homology over ℤ, computed by
[SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl). It has been validated on the 344 reference resolutions of
the preliminary implementation, among which those of the dual braid monoids of the 26
well-generated exceptional complex reflection groups: the integral homology of the corresponding complex braid groups is reproduced. The construction of the resolutions themselves is being
ported, with integer-coded simples and precomputed lcm tables.

```julia
using Gauss, SparseArrays

# Cellular chain complex of the real projective plane: d₀ = 0, d₁ = 2.
C = ChainComplex([spzeros(Int, 1, 1), sparse([1], [1], [2], 1, 1)])
homology(C)    # ℤ, ℤ/2ℤ, 0 in degrees 0 to 2
```

## Requirements

Julia 1.11 or later. Monoids are built from reflection groups through
[Chevie.jl](https://github.com/jmichel7/Chevie.jl), loaded as a package extension.
Elementary divisors of the large sparse differentials are computed by
[SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl).

## Development

This package is an overhauled version of the code available here: [Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl). The old version was developed by the author alone, while the new version is developed with the assistance of an AI model (Claude Opus 5.5). All code is reviewed, tested and validated by the author against reference computations.

## License

MIT, see [LICENSE](LICENSE).
