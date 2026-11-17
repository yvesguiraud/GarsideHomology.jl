# Complexes

```@meta
CurrentModule = GarsideHomology.Complexes
```

```@docs
Complexes
```

The exported names of `GarsideHomology.Complexes` are also exported by `GarsideHomology`. Some of its public names are also public in `GarsideHomology`: `rank`, `basis`, `MonoidAlgebra`, `GradedBasis`, `BasedFreeModuleElement` and `Differential`. The others are written with the name of the submodule, as in `GarsideHomology.Complexes.AbstractChainComplex`. With Chevie loaded, `basis` and `rank` can be written without qualification, since the extension adds methods to the functions of Chevie that have these names.

## Rings and monoid algebras

```@docs
AbstractMonoid
MonoidAlgebra
```

## Modules

```@docs
AbstractModule
BasedFreeModule
basis
```

## Labelled sparse matrices

```@docs
LabelledSparseMatrix
labels
extend_scalars
```

## Chain complexes

The docstrings of `AbstractChainComplex` include the methods that build the elements of a complex: `C("…")` from their written form, and, for a resolution, `P(x)` from a cell and `P(G(2), G(1, 2))` from germ elements. The cells themselves are given by [`Cell`](@ref GarsideHomology.Cell).

```@docs
AbstractChainComplex
ChainComplex
differential
is_square_zero
collapsible_pairs
```

## Graded bases and elements

```@docs
GradedBasis
BasedFreeModuleElement
Differential
```

## Display

```@docs
Notation
notation
set_notation!
```

## Finitely generated modules

```@docs
FinitelyGeneratedModule
FinitelyGeneratedAbelianGroup
torsion
FinitelyGeneratedGradedModule
FinitelyGeneratedGradedAbelianGroup
```

The function `homology` is documented, with all its methods, on the page [Coefficients and homology](coefficients.md).

## For developers

```@docs
label_type
identity_label
check_invariants
```
