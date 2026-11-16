# Germs

```@meta
CurrentModule = GarsideHomology
```

A germ gives the simples of a Garside monoid, numbered by integers. For simples `s` and `t` and an atom `a`, it gives:
- the product `st`, when it is simple;
- the right complement `s \ t`, that is, the simple `u` such that `su = t`, when `s ≼ t`;
- the right least common multiple (right lcm) `s ∨ a`.

Exported: `natoms`, `nsimples`, `simples`, `opposite`, `germ`, `dual_germ`. The other names are public only. With Chevie loaded, `word` and `elements` can be written without qualification, since the extension adds methods to the functions of Chevie that have these names.

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
