# Two orderings

```@meta
CurrentModule = GarsideHomology
```

The ordering of an lcm resolution ([`LCMOrdering`](@ref)) fixes its cells: at each node of the tree of the cells, a preorder on the lcms with the atoms decides which atoms give children (see [Background and conventions](../user_guide/background.md)). This page compares the two main orderings on the dual braid monoid B*(B₃), whose 9 atoms are the reflections of the Coxeter group B₃:
- [`ComplementDeglex`](@ref), the default, which gives the reduced resolution;
- [`TrivialPreorder`](@ref), with the order `1:9` of the atoms at the root, which gives the order resolution of Dehornoy and Lafont for this order.

The numbers and the trees below are computed by the package when this site is built.

## Numbers of cells

```@example orderings
using GarsideHomology, Chevie

G = dual_germ(coxgroup(:B, 3))
R = lcm_resolution(G)                                    # ComplementDeglex
D = lcm_resolution(G, LCMOrdering(TrivialPreorder()))    # Dehornoy–Lafont, order 1:9
[ncells(R, n) for n in 0:3], [ncells(D, n) for n in 0:3]
```

The reduced resolution has 38 cells, and the order resolution 40: one more cell in degree 2, and one more in degree 3. The reduced resolution is minimal, so that its numbers of cells do not depend on the choices of the construction.

## Two fragments of the trees

The function below prints the subtree of a cell, with the atom that created each cell. It reads the tree of the cells through its accessors ([`CellTree`](@ref)).

```@example orderings
function print_subtree(T, x, depth = 0)
    atom = x == 1 ? "" : "   (atom $(G(GarsideHomology.cell_atom(T, x))))"
    println("  "^depth, GarsideHomology.Cell(T, Int32(x)), atom)
    for c in GarsideHomology.cell_children(T, x)
        print_subtree(T, c, depth + 1)
    end
end

# The subtree of the cell of degree 1 whose letter is the atom a, for both orderings.
function compare_subtrees(a)
    for (name, P) in (("ComplementDeglex", R), ("TrivialPreorder", D))
        T = GarsideHomology.cell_tree(P)
        x = only(x for x in GarsideHomology.cell_range(T, 1) if GarsideHomology.bar_form(T, x) == [a])
        println(name, ":")
        print_subtree(T, x, 1)
    end
end
nothing # hide
```

The children of the cell `9` come from the atoms 1, 2 and 3, whose lcms with `9` are `9·48`, `9·8` and `9·4`:

```@example orderings
compare_subtrees(9)
```

- `ComplementDeglex` sorts these atoms by the right complements of their lcms, the shorter first: `4` (atom 3), `8` (atom 2), then `48` (atom 1). An atom is kept only if no atom placed before it divides its lcm; the atom 1 is therefore discarded, since `9 ∨ 3 = 9·4` divides `9 ∨ 1 = 9·48`.
- `TrivialPreorder` keeps the order of the root, 1, 2, 3: the atom 1 comes first and is kept, which gives the cell `9|48`, whose last letter is not an atom. The list of admissible atoms of a child is formed of the atoms of its elder siblings, so that the cell `9|4` now has the list `[1, 2]`, which gives the cell `9|4|8`.

The order of the children also matters when no atom is discarded. The children of the cell `4` come from the atoms 1 and 3, with the complements `6` and `3`:

```@example orderings
compare_subtrees(4)
```

Both orderings keep both atoms, but not in the same order. With `ComplementDeglex`, the atom 3 comes first, so that the second child, `4|6`, inherits the list `[3]` and has the child `4|6|3`. With `TrivialPreorder`, the atom 1 comes first, and the second child, `4|3`, inherits the list `[1]` and has the child `4|3|8`. The numbers of cells are the same, but the cells differ.
