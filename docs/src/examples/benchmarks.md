# Benchmarks

```@meta
CurrentModule = GarsideHomology
```

This page gives the times of the construction of the reduced lcm resolution of seven Garside monoids, from the germ to the differentials, and compares them with those of the preliminary implementation [Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl). The measurements were made on 18 November 2026 on a workstation with an Intel Core Ultra 9 285 processor (24 threads, 128 GB of memory, Windows), with Julia 1.13.4 started with all its threads, Chevie 0.1.15 and the commit 169582c of the package. They were made by the scripts of the folder `benchmark` of the repository, where the raw results are kept (`benchmark/results`).

## Protocol

The protocol is the following:
- each run is made in a fresh Julia session;
- a warm-up on a small monoid of the same type compiles the code first, so that the times below do not include the compilation;
- each time is the minimum of three runs, and each peak memory is the largest of the three; the three runs of a case differ by at most 5%;
- the construction of the reflection group by Chevie is not counted.

The peak memory is that of the whole Julia process, which takes about 0.8 GiB once Julia, Chevie and the package are loaded. The package builds the resolutions on a single thread: the number of threads plays no role here.

## Construction with the package

The package builds the resolution in three steps, which are timed separately:
1. the germ, by [`dual_germ`](@ref) for the dual braid monoids, and by [`germ`](@ref) for the classical braid monoid B⁺(E₈), whose germ is lazy (a [`CoxeterGerm`](@ref)), since the Coxeter group E₈ has 696729600 elements;
2. the tree of the generating cells ([`CellTree`](@ref));
3. the differentials ([`LCMProjectionData`](@ref) and [`lcm_resolution`](@ref)).

All the dual braid monoids are built with the Coxeter element given by the product of the generators of the group in the order of Chevie. The groups G₃₅, G₃₆ and G₃₇ are the Coxeter groups E₆, E₇ and E₈.

| Monoid | Atoms | Simples | Generating cells | Terms | Germ | Tree | Differentials | Total | Peak memory |
|---|---|---|---|---|---|---|---|---|---|
| B*(G₃₄) | 56 | 1584 | 12796 | 141792 | 0.21 s | < 0.01 s | 0.024 s | 0.23 s | 0.86 GiB |
| B*(G₃₅) | 36 | 833 | 4878 | 48334 | 0.092 s | < 0.01 s | < 0.01 s | 0.10 s | 0.79 GiB |
| B*(G₃₆) | 63 | 4160 | 39440 | 488741 | 0.032 s | < 0.01 s | 0.13 s | 0.17 s | 0.87 GiB |
| B*(G₃₇) | 120 | 25080 | 390720 | 6071842 | 0.41 s | 0.014 s | 3.25 s | 3.68 s | 1.75 GiB |
| B*(H₄) | 60 | 280 | 1080 | 7795 | 0.26 s | < 0.01 s | < 0.01 s | 0.26 s | 0.81 GiB |
| B*(D₈) | 56 | 9438 | 108976 | 1486032 | 0.12 s | < 0.01 s | 0.67 s | 0.79 s | 1.04 GiB |
| B⁺(E₈) | 8 | 696729600 | 256 | 909652 | < 0.01 s | < 0.01 s | 1 min 14 s | 1 min 14 s | 9.13 GiB |

B⁺(E₈) is the most expensive case. Since its germ is lazy, the products of simples are computed during the construction of the differentials, where the projection onto the generating cells visits 41 million chains of the Garside complex to give the 909652 terms of the differentials.

The number of terms depends on the numbering of the atoms, which is not always the same in `dual_germ` and in the monoids of Chevie: the table below, made with the monoids of Chevie, gives other numbers of terms for B*(G₃₄) to B*(G₃₇), with the same numbers of generating cells.

## Comparison with Gauss.jl

[Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl) is the preliminary implementation of the lcm resolutions, which the package replaces; its variant `Gauss_a` gives the reduced lcm resolution. It works on the monoids built by Chevie. For the comparison, the package is run on the same monoids, through [`germ`](@ref): it then builds exactly the same resolutions as Gauss.jl, with the same generating cells and the same terms.

The times of Gauss.jl are those of its two steps, the generating cells and the differentials (the functions `lcm_gen_data` and `filtered_lcm_differential`), and those of the package include the germ, the tree and the differentials. Neither includes the construction of the monoid by Chevie, which takes 4 min 1 s for B*(G₃₄) and less than 0.05 s for the other monoids; the peak memories include it. The runs of Gauss.jl were made case by case, the three runs of a case one after the other, while those of the package were made in random order over all the cases.

| Monoid | Generating cells | Terms | Gauss.jl | Package | Ratio | Peak memory, Gauss.jl | Peak memory, package |
|---|---|---|---|---|---|---|---|
| B*(G₃₄) | 12796 | 144406 | 42 min 43 s | 0.17 s | 14900 | 2.12 GiB | 2.43 GiB |
| B*(G₃₅) | 4878 | 50411 | 3.70 s | 0.065 s | 57 | 0.83 GiB | 0.84 GiB |
| B*(G₃₆) | 39440 | 539572 | 22.5 s | 0.23 s | 98 | 1.79 GiB | 0.89 GiB |
| B*(G₃₇) | 390720 | 7005479 | 7 min 11 s | 4.09 s | 105 | 24.3 GiB | 1.99 GiB |
| B*(H₄) | 1080 | 7795 | 1 min 55 s | 0.22 s | 523 | 0.88 GiB | 0.84 GiB |
| B*(D₈) | 108976 | 1486032 | 7 min 47 s | 0.84 s | 556 | 3.75 GiB | 1.05 GiB |
| B⁺(E₈) | 256 | 909652 | 19 min 20 s | 1 min 13 s | 15.9 | 41.7 GiB | 8.94 GiB |

The ratio is the time of Gauss.jl divided by that of the package. The computation of the homology from the differentials is done by [SparseElementaryDivisors.jl](https://github.com/yvesguiraud/SparseElementaryDivisors.jl), whose documentation gives times for the differentials of B*(G₃₇).
