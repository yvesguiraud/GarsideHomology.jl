# Lines of the snapshots, format version 1 (test/data/reference/FORMAT.md). The lines of
# the cells (C) come with the cell tree (J1.d, task d2), those of the differential (D)
# with the projection (task d3); the header and the writing of the files come with the
# resolution (task d4).

# A simple by its smallest word, the atoms separated by dots; the identity is "e".
_simple_string(G::AbstractGerm, s::Integer) = iszero(s) ? "e" : join(word(G, s), '.')

"""
    _cell_lines(T::CellTree, n; bar = false) -> Vector{String}

The lines `C <degree> <atom sequence> <bar form> <weight>` of the cells of degree `n` of
`T`, in the order of the cells (to be sorted by bytes for the fingerprints); with
`bar = true`, the lines `C <degree> <bar form> <weight>` of the fingerprints in bar form.
"""
function _cell_lines(T::CellTree, n::Integer; bar::Bool = false)
    G = T.germ
    words = Dict{Int32, String}()
    simple(s) = get!(() -> _simple_string(G, s), words, Int32(s))
    lines = String[]
    for x in cell_range(T, n)
        b = n == 0 ? "[]" : "[" * join((simple(s) for s in bar_form(T, x)), '|') * "]"
        w = simple(cell_weight(T, x))
        if bar
            push!(lines, "C $n $b $w")
        else
            a = n == 0 ? "-" : join(atom_sequence(T, x), '.')
            push!(lines, "C $n $a $b $w")
        end
    end
    return lines
end

"""
    _differential_lines(P::LCMProjectionData, n; bar = false) -> Vector{String}

The lines `D <n> <source> <coef> <σ> <target> <τ>` of the differential `dₙ`, from the
cells of degree `n + 1` to those of degree `n`: one line for each term `coef · σ [target] τ`
of the boundary of a source, the cells named by their atom sequences (`-` for the cell of
degree 0), or by their bar forms with `bar = true`; the coefficient with its sign (`+1`,
`-2`), the simples by their smallest words. In the order of the cells (to be sorted by
bytes for the fingerprints).
"""
function _differential_lines(P::LCMProjectionData, n::Integer; bar::Bool = false)
    T = P.tree
    G = T.germ
    words = Dict{Int32, String}()
    simple(s) = get!(() -> _simple_string(G, s), words, Int32(s))
    name(x) = bar ? "[" * join((simple(s) for s in bar_form(T, x)), '|') * "]" :
              (cell_degree(T, x) == 0 ? "-" : join(atom_sequence(T, x), '.'))
    lines = String[]
    for x in cell_range(T, n + 1)
        sx = name(x)
        for t in lcm_boundary(P, x)
            c = t.coef > 0 ? "+$(t.coef)" : "-$(-t.coef)"
            push!(lines, "D $n $sx $c $(simple(t.left)) $(name(t.cell)) $(simple(t.right))")
        end
    end
    return lines
end
