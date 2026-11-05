# Complexes

```@meta
CurrentModule = GarsideHomology.Complexes
```

```@docs
Complexes
```

The exported names of `GarsideHomology.Complexes` are also exported by `GarsideHomology`. Among its public names, `rank`, `basis`, `MonoidAlgebra`, `GradedBasis`, `BasedFreeModuleElement` and `Differential` are also public in `GarsideHomology`; the others are written `GarsideHomology.Complexes.AbstractChainComplex`, for instance. With Chevie loaded, `basis` and `rank` are written without qualification: the extension adds methods to the functions of Chevie that have these names.

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

The docstrings of `AbstractChainComplex` include the methods that build the elements of a complex: `C("…")` from their written form, and, for a resolution, `P(x)` from a cell and `P(G(2), G(1, 2))` from germ elements.

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

The homology, `homology`, is documented on the page [Coefficients and homology](coefficients.md), with all its methods.

## For developers

```@docs
label_type
identity_label
check_invariants
```
