# GarsideHomology.jl

Small free resolutions of Garside monoids and their homology, in Julia.

!!! note "Status"
    Early development, not registered yet: the interface may still change.

## Introduction

### Garside monoids and braid groups

A *Garside monoid* is a cancellative monoid in which divisibility is a lattice order, with a *Garside element* Δ whose divisors, the *simples*, are finite in number and generate the monoid. Its enveloping group is a group of fractions, called a *Garside group*, and it has the same homology as the monoid. Through their classical and dual braid monoids, the braid groups of the finite real reflection groups and of the well-generated complex reflection groups are Garside groups. The notions and the conventions of the package are recalled on the page [Background and conventions](user_guide/background.md).

### Computation method

GarsideHomology.jl computes the *reduced lcm resolution* of a Garside monoid, a variant of the resolution of Dehornoy and Lafont with fewer cells. It is a free resolution of the monoid algebra by free bimodules, whose cells are chains of simples, and it can be checked by a certificate: homogeneity of the boundaries, d² = 0, and Euler characteristic of the cells of each weight. The homology of the group with coefficients in a module of rank one (for now, the trivial module ℤ and the sign module ℤε, in which the braid reflections act by −1) is then computed from the elementary divisors of the differentials, by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl).

### Motivation

The first target is the homology of the braid groups of the real and well-generated complex reflection groups. For the 26 well-generated exceptional groups, the package computes the homology with ℤ and with ℤε, except for B(G₃₇) with ℤε, where the coefficients of the elimination outgrow the integer types used for now. The values agree with the tables of Garnier, with one exception noted on the page [Exceptional braid groups](examples/exceptional.md).

## Quick start

### Installation

The package is not registered yet:

```julia
using Pkg
Pkg.add(url = "https://github.com/yvesguiraud/GarsideHomology.jl")
Pkg.add("Chevie")   # optional: monoids built from reflection groups
```

It requires Julia 1.11 or later. The monoids are built from the reflection groups through [Chevie.jl](https://github.com/jmichel7/Chevie.jl), loaded as a package extension.

### A first example

```julia
using GarsideHomology, Chevie

G = dual_germ(crg(4))          # the germ of the dual braid monoid B*(G₄)
R = lcm_resolution(G)          # its reduced lcm resolution
certify(R)                     # correct, minimal, linear
homology(R, ℤ)                 # ℤ, ℤ, 0 in degrees 0 to 2
homology(R, ℤε)                # ℤ/2ℤ, ℤ/3ℤ, 0
```

The braid group of G₄ is the braid group on three strands. Start Julia with several threads to benefit from the parallel steps: `julia -t auto`. The [User guide](user_guide/getting_started.md) explains how to use the package, from this example on, and the [Reference manual](reference_manual/overview.md) gathers the docstrings, layer by layer.

## About

### Development

This package is developed with the help of an AI assistant. In every version, the part of the code that is produced by the AI assistant is reviewed by the author and validated by tests against reference computations.

GarsideHomology.jl is an overhauled version of the code available at [Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl), which was developed by the author alone; the new version is developed with Claude, currently with the Claude Opus 5.5 model.

### Acknowledgements

The author thanks Jean Michel for having introduced him to Julia a few years ago, and for his regular and valuable advice on the language. His package [Chevie.jl](https://github.com/jmichel7/Chevie.jl) is particularly useful for generating concrete test examples.

### Licence

MIT, see the file [LICENSE](https://github.com/yvesguiraud/GarsideHomology.jl/blob/main/LICENSE) of the repository.

### Contact

Feedback is most welcome, by [email](mailto:yves.guiraud@imj-prg.fr) to the author, or as an [issue](https://github.com/yvesguiraud/GarsideHomology.jl/issues) on the repository.
