# Orderings

```@meta
CurrentModule = GarsideHomology
```

An ordering of the lcms chooses the resolution: the reduced one with `ComplementDeglex()` (the default), that of Dehornoy and Lafont with `TrivialPreorder()`.

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

The vocabulary for writing a new preorder.

```@docs
AbstractLCMPreorder
precedes
extends_divisibility
```
