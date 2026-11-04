# Germs

```@meta
CurrentModule = GarsideHomology
```

A germ gives the simples of a Garside monoid, numbered by integers, with their products and their left complements.

Exported: `natoms`, `nsimples`, `simples`, `opposite`, `germ`, `dual_germ`. The other names are public only. With Chevie loaded, `word` and `elements` are written without qualification: the extension adds methods to the functions of Chevie that have these names.

## Interface of a germ

```@docs
AbstractGerm
natoms
nsimples
simples
elements
Base.one
Base.isone
atoms
product
ldiv
is_left_divisor
left_descents
lcm_atom
opposite
word
isless_lex
isless_deglex
```

## Elements

```@docs
GermElement
```

## Construction

```@docs
germ
dual_germ
TableGerm
CoxeterGerm
hurwitz_orbit
```

## Checks

```@docs
verify
count_multichains
GermFingerprint
fingerprint
germ_isomorphism
```

## For developers

```@docs
StaticBitSet
tabulate
nnumbered
```
