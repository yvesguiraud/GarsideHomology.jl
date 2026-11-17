# GarsideHomology.jl

[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yvesguiraud.github.io/GarsideHomology.jl/dev/)
[![CI](https://github.com/yvesguiraud/GarsideHomology.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yvesguiraud/GarsideHomology.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/yvesguiraud/GarsideHomology.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yvesguiraud/GarsideHomology.jl)

Small free resolutions of Garside monoids and their homology, in Julia.

GarsideHomology.jl computes the *reduced lcm resolution* of a Garside monoid, a variant of the resolution of Dehornoy and Lafont with fewer cells. It is a free resolution of the monoid algebra by free bimodules, and it can be checked by a certificate. The homology of the group with coefficients in a module of rank one (for now, the trivial module ℤ and the sign module ℤε, in which the braid reflections act by −1) is then computed from the elementary divisors of the differentials, by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl).

The first target is the homology of the braid groups of the real and well-generated complex reflection groups. For the 26 well-generated exceptional groups, the package computes the homology with ℤ and with ℤε, except for B(G₃₇) with ℤε, where the coefficients of the elimination outgrow the integer types used for now. The values agree with the tables of Garnier, with one exception noted on the page [Exceptional braid groups](https://yvesguiraud.github.io/GarsideHomology.jl/dev/examples/exceptional/) of the documentation.

> **Status:** early development, not registered yet; the interface may still change.

The [documentation](https://yvesguiraud.github.io/GarsideHomology.jl/dev/) contains a user guide, examples, the reference manual and a bibliography, as well as the development history and the acknowledgements.

## Installation

The package is not registered yet:

```julia
using Pkg
Pkg.add(url = "https://github.com/yvesguiraud/GarsideHomology.jl")
Pkg.add("Chevie")   # optional: monoids built from reflection groups
```

It requires Julia 1.11 or later. The monoids are built from the reflection groups through [Chevie.jl](https://github.com/jmichel7/Chevie.jl), loaded as a package extension.

## A first example

```julia
using GarsideHomology, Chevie

G = dual_germ(crg(4))          # the germ of the dual braid monoid B*(G₄)
R = lcm_resolution(G)          # its reduced lcm resolution
certify(R)                     # correct, minimal, linear
homology(R, ℤ)                 # ℤ, ℤ, 0 in degrees 0 to 2
homology(R, ℤε)                # ℤ/2ℤ, ℤ/3ℤ, 0
```

Start Julia with several threads to benefit from the parallel steps: `julia -t auto`. The page [Working with other packages](https://yvesguiraud.github.io/GarsideHomology.jl/dev/user_guide/interoperability/) of the user guide explains how to type `ℤ` and `ℤε`, which names are shared with Oscar, Hecke and Chevie, and the current limits of the package.
