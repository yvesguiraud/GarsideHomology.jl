# Cells

```@meta
CurrentModule = GarsideHomology
```

The cells of an lcm resolution are the nodes of a tree, numbered by integers. A [`Cell`](@ref) holds a number together with its tree.

Exported: `ncells`, `cell_weight`, `cell_degree`, `atom_sequence`. The other names are public only.

## Cells

```@docs
Cell
cell_weight(::Cell)
cell_degree(::Cell)
atom_sequence(::Cell)
```

## Tree of the cells

```@docs
CellTree
ncells
```

## For developers

### Accessors of the tree

These functions read the tree by the numbers of the cells.

```@docs
cell_range
cell_weight(::CellTree, ::Integer)
cell_degree(::CellTree, ::Integer)
cell_letter
cell_atom
bar_form
atom_sequence(::CellTree, ::Integer)
cell_children
admissible_atoms
```

### Projection

The projection computes the boundaries of the cells. Its names are provisional.

```@docs
ChainTree
BimoduleTerm
LCMProjectionData
lcm_splitting
lcm_projection!
lcm_boundary
```
