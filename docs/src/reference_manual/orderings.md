# Orderings

```@meta
CurrentModule = GarsideHomology
```

An ordering of the lcms chooses the resolution: `ComplementDeglex()`, the default, gives the reduced resolution, and `TrivialPreorder()` gives that of Dehornoy and Lafont.

Exported: `LCMOrdering`, `ComplementDeglex`, `TrivialPreorder`, `FunctionPreorder`. Public only: the names of the section For developers.

## Orderings

```@docs
LCMOrdering
```

## Preorders

```@docs
ComplementDeglex
TrivialPreorder
FunctionPreorder
```

## For developers

These names serve to write a new preorder.

```@docs
AbstractLCMPreorder
precedes
extends_divisibility
```
