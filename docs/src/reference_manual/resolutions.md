# Resolutions and certificates

```@meta
CurrentModule = GarsideHomology
```

Exported: `lcm_resolution`, `certify`. The other names are public only.

## Resolutions

```@docs
lcm_resolution
LCMResolution
CellBasis
cell_tree
```

A cell of a resolution `P` is given by its written form, `Cell(P, "2|12")`, or by its letters, `Cell(P, G(2), G(1, 2))` (see [`Cell`](@ref)). Calling the resolution gives the elements of its modules: `P(x)` is the generator `[x]` of a cell `x`, `P("2|12")` and `P(G(2), G(1, 2))` are the generator `[2|12]`, and `P("2×1[2]1 - [1]")` reads an element from its written form. These methods, and the differential, are documented with [`AbstractChainComplex`](@ref), on the page [Complexes](complexes.md).

## Certificates

```@docs
certify
ResolutionCertificate
is_correct
is_minimal
is_linear
```

The check of d² = 0 alone, [`is_square_zero`](@ref), is documented on the page [Complexes](complexes.md).

## For developers

```@docs
write_snapshot
```
