# Coefficients and homology

```@meta
CurrentModule = GarsideHomology
```

Exported: `TrivialModule`, `SignModule`, `ℤ`, `ℤε`, `ℤϵ`, `homology`. Public only: `AbstractMonoidModule`, `tensor_product`, `⊗`, and `ZZ`, the ring of integers of AbstractAlgebra, which is the same object as `ℤ`.

Some of these names are also those of other packages:
- `ZZ` is not exported, because Nemo, Hecke and Oscar export a `ZZ` of their own;
- with Hecke or Oscar loaded, `homology` and `differential` must be written `GarsideHomology.homology` and `GarsideHomology.differential`, since these packages also export functions with these names;
- with Chevie loaded, `⊗` can be written without qualification, since the extension adds methods to the `⊗` of Chevie.

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

This section gathers all the methods of `homology`: for a chain complex over ℤ, for a resolution with coefficients, and for a germ.

```@docs
homology
```
