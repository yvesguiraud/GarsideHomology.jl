# Exceptional braid groups

```@meta
CurrentModule = GarsideHomology
```

The 26 well-generated exceptional groups, in the order of the tables of Garnier. The column *Class* gives the real type of W, or the representative of its isodiscriminantality class: the braid groups of a class are isomorphic, by an isomorphism that sends braid reflections to braid reflections, hence have the same homology for both modules.

Trivial coefficients:

```@eval
Main.ExamplesTables.exceptional_homology("trivial")
```

Sign module ℤε:

```@eval
Main.ExamplesTables.exceptional_homology("sign")
```

† The value with ℤε for G₃₇ is not computed by the package for now: the coefficients of the elimination outgrow `Int128`, the largest integer type used in these computations. The value given was obtained by computing the elementary divisors of the differentials of the package with Magma.

The eight other exceptional groups, G₇, G₁₁, G₁₂, G₁₃, G₁₅, G₁₉, G₂₂ and G₃₁, are not well-generated: they have no dual braid monoid, and the package does not compute their homology for now. Their homology with ℤ and ℤε is known: see Callegaro and Marin, Marin (*Homology computations for complex braid groups II*) and Garnier (tables 3 and 4), who treats G₃₁ through a Garside category ([Bibliography](../bibliography.md)).

Numbers of cells of the reduced lcm resolution:

```@eval
Main.ExamplesTables.exceptional_cells()
```

For five of these groups, table 1 of Garnier gives the numbers of cells of Dehornoy–Lafont complexes, for the order of the atoms of Marin and for an order chosen by Garnier to reduce the number of cells (with the correction of his erratum):

```@eval
Main.ExamplesTables.exceptional_dehornoy_lafont()
```
