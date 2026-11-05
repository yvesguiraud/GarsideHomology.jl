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

The elements of a resolution `P` are built by `P(x)` for a cell `x`, `P(G(2), G(1, 2))` from germ elements, and `P("2|12")` from their written form: these methods are documented with [`AbstractChainComplex`](@ref), on the page [Complexes](complexes.md), as is the differential.

## Certificates

```@docs
certify
ResolutionCertificate
is_correct
is_minimal
is_linear
```

The check of ∂² = 0 alone, [`is_square_zero`](@ref), is on the page [Complexes](complexes.md).

## For developers

```@docs
write_snapshot
```
