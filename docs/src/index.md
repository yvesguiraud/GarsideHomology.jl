# GarsideHomology.jl

Small free resolutions of Garside monoids and their homology, in Julia.

!!! note "Status"
    Early development, not registered yet: the interface may still change.

## Introduction

### Garside monoids and braid groups

A *Garside monoid* is a cancellative monoid in which divisibility is a lattice order, and which has a distinguished element, called the *Garside element*, whose divisors, called the *simples*, are finite in number and generate the monoid. A *Garside group* is the enveloping group of a Garside monoid; it is a group of fractions, and it has the same homology as the monoid. The seminal example of a Garside monoid and of a Garside group are the usual braid monoid and the usual braid group; more generally, the classical and dual braid monoids of the finite real reflection groups, and the dual braid monoids of the well-generated complex reflection groups, are Garside monoids, while the corresponding braid groups are Garside groups. The notions and the conventions of the package are recalled on the page [Background and conventions](user_guide/background.md).

### Computation method

GarsideHomology.jl computes *resolutions by least common multiples*, or *lcm resolutions*, of Garside monoids: this is a family of resolutions of the monoid algebra by free bimodules, whose generators are chains of simples corresponding to lcms of some sequences of atoms. An lcm resolution is parameterised by an order on the atoms and their lcms; among the possible choices of orders, a trivial one recovers the order resolution of Dehornoy and Lafont [DL99, §4], while another choice gives a *reduced lcm resolution* with fewer generating cells than the order resolutions.

The construction of the resolution can be checked by certificates: differential squaring to zero, homogeneity of the boundaries, and Euler characteristic of the generating cells. The homology of the group with coefficients in a module of rank one (for now, the trivial module ``ℤ`` and the sign module ``ℤ_ε``, in which the braid reflections act by ``−1``) is then computed from the elementary divisors of the differentials, by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl).

### Motivation

The first target is the homology of the braid groups of the real and well-generated complex reflection groups (see Examples). For the 26 well-generated exceptional groups, the package computes the homology with ``ℤ`` and with ``ℤ_ε``, except for ``B(G₃₇)`` with ``ℤ_ε`` where the coefficients of the elimination outgrow the integer types used for now; the values agree with the tables of [Gar24, §3], with one exception noted on the page [Exceptional braid groups](examples/exceptional.md).

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

W = crg(4)                     # the complex reflection group G₄ (constructed by Chevie)
G = dual_germ(W)               # the germ of the dual braid monoid B^*(G₄)
P = lcm_resolution(G)          # its reduced lcm resolution
certify(P)                     # correct, minimal, linear
homology(P, ℤ)                 # ℤ, ℤ, 0 in degrees 0 to 2
homology(P, ℤε)                # ℤ/2ℤ, ℤ/3ℤ, 0 in degrees 0 to 2
```

Start Julia with several threads to benefit from the parallel steps: `julia -t auto`. The [User guide](user_guide/getting_started.md) explains how to use the package, from this example on, and the [Reference manual](reference_manual/overview.md) gathers the docstrings, layer by layer.

## About

### Development

GarsideHomology.jl is an overhauled version of the code available at [Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl), which was developed by the author alone. The present package is developed with the help of an AI assistant (Claude, currently with the Opus 5.5 model). In every version, the part of the code that is produced by the AI assistant is reviewed by the author and validated by tests against reference computations.

### Acknowledgements

The author thanks Jean Michel for having introduced him to Julia a few years ago, and for his regular and valuable advice on the language; his package [Chevie.jl](https://github.com/jmichel7/Chevie.jl) is particularly useful for generating concrete test examples. He also thanks Owen Garnier for many fruitful conversations about the Dehornoy–Lafont's resolution, and in particular on the possible choices of orders of atoms that led to the discovery of the reduced lcm resolution.

### Licence

MIT, see the file [LICENSE](https://github.com/yvesguiraud/GarsideHomology.jl/blob/main/LICENSE) of the repository.

### Contact

Feedback is most welcome, by [email](mailto:yves.guiraud@imj-prg.fr) to the author, or as an [issue](https://github.com/yvesguiraud/GarsideHomology.jl/issues) on the repository.
