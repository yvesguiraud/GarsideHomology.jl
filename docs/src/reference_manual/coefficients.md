# Coefficients and homology

```@meta
CurrentModule = GarsideHomology
```

Exported: `TrivialModule`, `SignModule`, `ℤ`, `ℤε`, `ℤϵ`, `homology`. Public only: `AbstractMonoidModule`, `tensor_product`, `⊗`, and `ZZ`, the ring of integers of AbstractAlgebra, the same object as `ℤ`. `ZZ` is not exported because Nemo, Hecke and Oscar export a `ZZ` of their own; with Hecke or Oscar loaded, `homology` and `differential` are to be written `GarsideHomology.homology` and `GarsideHomology.differential`, since they also export functions with these names. With Chevie loaded, `⊗` is written without qualification: the extension adds methods to the `⊗` of Chevie.

## Modules of coefficients

```@docs
AbstractMonoidModule
TrivialModule
SignModule
ℤε
ℤϵ
```

## Rings

```@docs
ℤ
```

## Tensor products

```@docs
tensor_product
⊗
```

## Homology

All the methods of `homology` are here: of a chain complex over ℤ, of a resolution with coefficients, and of a germ.

```@docs
homology
```
