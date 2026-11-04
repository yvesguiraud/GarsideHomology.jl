# Cells

```@meta
CurrentModule = GarsideHomology
```

The cells of an lcm resolution are the nodes of a tree, numbered by integers; a [`Cell`](@ref) carries a number together with its tree.

Exported: `ncells`, `cell_weight`, `cell_degree`, `atom_sequence`. The other names are public only.

## Cells

```@docs
Cell
cell_weight
cell_degree
atom_sequence
```

## Tree of the cells

```@docs
CellTree
ncells
```

## For developers

The accessors of the tree, by cell number:

```@docs
cell_range
cell_letter
cell_atom
cell_children
admissible_atoms
bar_form
```

The projection, which computes the boundaries of the cells. Its names are provisional: they will follow the theory of the article in preparation.

```@docs
ChainTree
BimoduleTerm
LCMProjectionData
lcm_splitting
lcm_projection!
lcm_boundary
```
