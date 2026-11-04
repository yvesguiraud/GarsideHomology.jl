# Spherical Artin groups

```@meta
CurrentModule = GarsideHomology
```

The braid groups of the finite Coxeter groups, computed from the classical braid monoid B⁺(W) and from the dual braid monoid B*(W) wherever both were computed, with equal results. The types H₃, F₄, H₄, E₆, E₇ and E₈ are also G₂₃, G₂₈, G₃₀, G₃₅, G₃₆ and G₃₇ of the page [Exceptional braid groups](exceptional.md).

Trivial coefficients:

```@eval
Main.ExamplesTables.spherical_homology("trivial")
```

Sign module ℤε:

```@eval
Main.ExamplesTables.spherical_homology("sign")
```

Numbers of cells of the reduced lcm resolution of B*(W); that of B⁺(W), which is minimal too, has `binomial(n, k)` cells in degree `k`.

```@eval
Main.ExamplesTables.spherical_cells()
```

```@eval
Main.ExamplesTables.spherical_without_cells()
```

The dihedral types, from B*(W):

```@eval
Main.ExamplesTables.spherical_dihedral()
```
