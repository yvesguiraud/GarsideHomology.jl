# Examples

```@meta
CurrentModule = GarsideHomology
```

These pages gather what the package computes for now on the braid groups of reflection groups: the homology with trivial coefficients ℤ and with the sign module ℤε, in which the braid reflections act by −1, and the numbers of cells of the reduced lcm resolutions. They will be completed as the computations progress.

## How the values are obtained

For each group, the reduced lcm resolution of a Garside monoid of the braid group, then its certificate and its homology with ℤ and ℤε, the elementary divisors being computed by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl). The monoid is the dual braid monoid B*(W), for the Coxeter element chosen by default, and also the classical braid monoid B⁺(W) for the real groups.

The values are those of the files `test/data/expected/*.toml` of the [repository](https://github.com/yvesguiraud/GarsideHomology.jl). The long tests of the package recompute them (environment variable `GARSIDE_HOMOLOGY_LONG_TESTS`), and the tables of these pages are generated from these files at each build of the site. The groups are written torsion first, as in the tables of Garnier, and `?` marks a value not computed yet.

The numbers of cells are those of the reduced lcm resolution of B*(W) in each degree. In all these cases, this resolution is minimal, so that these numbers do not depend on the choices of the construction.

## The pages

- [Exceptional braid groups](exceptional.md): the 26 well-generated exceptional groups.
- [Spherical Artin groups](spherical.md): the braid groups of the finite Coxeter groups.
- [Imprimitive braid groups](imprimitive.md): the groups G(m, m, n) and G(m, 1, n).

## Not computed yet

With ℤε, the coefficients of the elimination outgrow `Int128`, the largest integer type used in these computations, for G₃₇ (see [Exceptional braid groups](exceptional.md)) and for:

```@eval
Main.ExamplesTables.imprimitive_unknown()
```

Not computed yet either: the eight exceptional groups that are not well-generated, the imprimitive groups beyond the ranges of [Imprimitive braid groups](imprimitive.md), and the homology with coefficients in the Laurent polynomials k[t, t⁻¹], in which the braid reflections act by t.
