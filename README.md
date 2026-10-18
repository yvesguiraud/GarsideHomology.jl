# Gauss.jl

[![CI](https://github.com/yvesguiraud/Gauss.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yvesguiraud/Gauss.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/yvesguiraud/Gauss.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yvesguiraud/Gauss.jl)

Small free resolutions of Garside monoids and their homology, in Julia.

Gauss.jl computes the *reduced lcm resolution* of a Garside monoid, a variant of the
Dehornoy–Lafont resolution with fewer cells, and the homology of the corresponding
groups with various coefficients. The first target is the homology of real and well-generated complex braid groups.

## Status

Early development, not registered. This version is a package skeleton: the computation
code is being ported from a preliminary implementation and rewritten with integer-coded
simples and precomputed lcm tables.

## Requirements

Julia 1.11 or later. Monoids are built from reflection groups through
[Chevie.jl](https://github.com/jmichel7/Chevie.jl), loaded as a package extension.
Elementary divisors of the large sparse differentials are computed by
[SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl).

## Development

This package is an overhauled version of the code available here: [Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl). The old version was developed by the author alone, while the new version is developed with the assistance of an AI model (Claude Opus 5.5).

## License

MIT, see [LICENSE](LICENSE).
