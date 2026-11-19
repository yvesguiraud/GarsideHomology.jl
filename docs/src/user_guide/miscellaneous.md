# Miscellaneous

```@meta
CurrentModule = GarsideHomology
```

## Display settings

The package provides settings to adjust the notations for atoms, their products, the elements of free modules, the differential, etc. The settings are stored using the [`Notation`](@ref) type: see there for the list of possible settings. The examples continue those of the page [Germs and resolutions](resolutions.md), in ``B⁺(A₂)``.

The settings can be changed for the whole session using [`set_notation!`](@ref), and `set_notation!()` goes back to the defaults:

```julia
2 * d(P("2"))                           # 2×2[] - 2×[]2
set_notation!(atom_names = ["a", "b"])
2 * d(P("2"))                           # 2b[] - 2[]b
P("2b[] - 2[]b") == 2 * d(P("2"))       # true: the reading follows the settings
set_notation!()
```

The settings can also be changed locally, using the `:notation` of an `IOContext`, taking precedence over the global configuration:

```julia
show(
    IOContext(stdout, :notation => Notation(coefficient_style = :parentheses)),
    2 * P("2")
)               # 2([2])
show(
    IOContext(stdout, :notation => Notation(empty_cell_in_module = "[e]")),
    d(P("1"))
)               # 1[e] - [e]1
```

## Names shared with other packages

Some names of GarsideHomology are also exported by other packages for computational algebra such as Nemo, Hecke and Oscar: `ZZ`, `homology`, `differential`, `tensor_product`, `⊗`, `elements`. To avoid conflicts about a name `xxx`, one can qualify it as in `GarsideHomology.xxx`, or import it explicitly to take precedence over other packages as in `using GarsideHomology: xxx`. Note that `ZZ`, `tensor_product` and `⊗` are public but not exported by GarsideHomology: they are always written `GarsideHomology.xxx` or imported, as in `using GarsideHomology: ZZ, ⊗` (the alias `ℤ` of `ZZ` is exported).

Some names are shared with Chevie but the extension GarsideHomologyChevieExt adds methods to Chevie functions to avoid conflicts: thus, `basis`, `word`, `rank`, `elements` and `⊗` can be used freely in that case. Without Chevie, `basis`, `word`, `rank` and `elements` are public but not exported: write `using GarsideHomology: basis, word, rank, elements`.

## Limits

- Chevie represents the reflection groups by permutations of their roots, with 16-bit integers by default, and the extension of GarsideHomology stores the images of the roots with 16-bit integers too: the groups with more than 32 767 roots, such as ``G(29, 29, 7)`` and its 35 322 roots, are out of reach for now.
- A tabulated germ ([`TableGerm`](@ref)) takes about 16 bytes per simple and per atom (49 MB for ``B^*(E₈)``). Above 5·10⁵ simples, the function `germ` gives a lazy germ for a classical braid monoid, but `dual_germ` throws an `ArgumentError`.
- A lazy germ ([`CoxeterGerm`](@ref)) numbers its simples as it meets them: it must not be queried from several threads at once.

## Snapshots

The function [`write_snapshot`](@ref)`(path, P)` writes a resolution in a compressed text format (`path` ending in `.txt.gz`), with its generating cells and its differential, and by default its fingerprints beside it (`.sha.csv`). This is the format of the reference resolutions of the tests; it is described in the file `test/data/reference/FORMAT.md` of the [repository](https://github.com/yvesguiraud/GarsideHomology.jl).