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

# ------------------------------------------------------------------------------------
# Writing (task d4)

_coefficient_string(c::Integer) = c > 0 ? "+$c" : "-$(-c)"

# The groups of lines of the snapshot of the lcm resolution C, in the order of the file:
# (section, n, lines in atom sequences, lines in bar forms), each sorted by bytes.
function _snapshot_groups(C::ChainComplex, T::CellTree)
    G = T.germ
    words = Dict{Int32, String}()
    simple(s) = get!(() -> _simple_string(G, s), words, Int32(s))
    atom_name(x) = cell_degree(T, x) == 0 ? "-" : join(atom_sequence(T, x), '.')
    bar_name(x) = "[" * join((simple(s) for s in bar_form(T, x)), '|') * "]"
    groups = NamedTuple{(:section, :n, :atoms, :bar), Tuple{String, Int, Vector{String}, Vector{String}}}[]
    for n in range(C)
        push!(groups, (section = "C", n = n, atoms = sort!(_cell_lines(T, n)),
                       bar = sort!(_cell_lines(T, n; bar = true))))
    end
    for n in first(range(C)):(last(range(C)) - 1)
        d = differential(C, n)
        source, target = cells(C[n + 1]), cells(C[n])
        rows, labs, vals = Complexes.SparseArrays.rowvals(d), labels(d), Complexes.SparseArrays.nonzeros(d)
        a = String[]
        b = String[]
        for j in 1:ncols(d)
            x = source[j]
            sa, sb = atom_name(x), bar_name(x)
            for k in Complexes.SparseArrays.nzrange(d, j)
                y = target[rows[k]]
                σ, τ = simple(labs[k][1]), simple(labs[k][2])
                c = _coefficient_string(vals[k])
                push!(a, "D $n $sa $c $σ $(atom_name(y)) $τ")
                push!(b, "D $n $sb $c $σ $(bar_name(y)) $τ")
            end
        end
        push!(groups, (section = "D", n = n, atoms = sort!(a), bar = sort!(b)))
    end
    return groups
end

# The value of the key empreinte_germe: computed from the germ (:compute), given (a
# GermFingerprint or a string), or omitted (nothing).
_fingerprint_string(G::AbstractGerm, f::Symbol) =
    f === :compute ? (applicable(simples, G) ? sprint(show, fingerprint(G)) : "non calculée (germe paresseux)") :
    throw(ArgumentError("write_snapshot: germ_fingerprint = :$f, expected :compute, nothing, a GermFingerprint or a string"))
_fingerprint_string(G::AbstractGerm, ::Nothing) = "non calculée"
_fingerprint_string(G::AbstractGerm, f::GermFingerprint) = sprint(show, f)
_fingerprint_string(G::AbstractGerm, f::AbstractString) = String(f)

function _snapshot_header(T::CellTree, case::AbstractString, extra, groups, germ_fingerprint)
    G = T.germ
    header = Pair{String, String}["format" => "cliche-gauss 1"]
    isempty(case) || push!(header, "cas" => String(case))
    push!(header, "ordre" => sprint(show, T.ordering))
    push!(header, "cote" => "ppcm à droite ; terme coef·σ·y·τ, σ à gauche, τ à droite ; d_n : C_{n+1} -> C_n")
    push!(header, "n_atomes" => string(natoms(G)))
    push!(header, "empreinte_germe" => _fingerprint_string(G, germ_fingerprint))
    for (k, v) in extra
        key, value = String(k), String(v)
        (any(c -> c == '\n' || c == '\r', key * value) || occursin(':', key)) &&
            throw(ArgumentError("write_snapshot: invalid header line $key: $value"))
        push!(header, key => value)
    end
    push!(header, "gauss" => string(pkgversion(@__MODULE__)))
    push!(header, "julia" => string(VERSION))
    push!(header, "date" => Libc.strftime("%Y-%m-%dT%H:%M:%S", time()))
    push!(header, "cellules" => string(sum(length(g.atoms) for g in groups if g.section == "C")))
    push!(header, "termes" => string(sum((length(g.atoms) for g in groups if g.section == "D"); init = 0)))
    return header
end

function _sha_lines(lines::Vector{String})
    ctx = SHA.SHA256_CTX()
    for l in lines
        SHA.update!(ctx, codeunits(l))
        SHA.update!(ctx, UInt8[0x0a])
    end
    return bytes2hex(SHA.digest!(ctx))
end

"""
    write_snapshot(path, C::ChainComplex; case = "", extra = Pair{String, String}[],
                   fingerprints = true, germ_fingerprint = :compute)

Write the lcm resolution `C` (over a monoid algebra of bimodules, with bases of cells:
see [`lcm_resolution`](@ref)) as a snapshot of format version 1, gzip-compressed, at
`path`, which must end with `.txt.gz`: the header, then the lines of the cells by
increasing degree and those of the differential by increasing `n`, sorted by bytes within
each group (test/data/reference/FORMAT.md). With `fingerprints = true`, also write the
fingerprints at the same place, with the extension `.sha.csv` (FORMAT.md, section 5).

The header gives the format, `case` (if not empty), the ordering, the side convention,
the number of atoms and the fingerprint of the germ ([`fingerprint`](@ref)), then the
pairs `key => value` of `extra` (what only the caller knows: label, construction,
machine, commits), then the versions of Gauss and Julia, the date and the numbers of
cells and terms (FORMAT.md, section 7). It does not enter the fingerprints.

The fingerprint of the germ can cost most of the writing for a large tabulated germ (such
as B⁺(A₈), whose 362 880 simples are tabulated): `germ_fingerprint` takes it
already computed (a `GermFingerprint` or a string), for instance once for the snapshots
of several orderings of the same germ, or `nothing` to omit it; by default (`:compute`),
it is computed if the germ lists its simples, and omitted for a lazy germ.
"""
function write_snapshot(path::AbstractString, C::ChainComplex; case::AbstractString = "",
                        extra = Pair{String, String}[], fingerprints::Bool = true,
                        germ_fingerprint = :compute)
    endswith(path, ".txt.gz") || throw(ArgumentError("write_snapshot: the path $path does not end with .txt.gz"))
    C.ring isa MonoidAlgebra && C.ring.side === :both ||
        throw(ArgumentError("write_snapshot: the complex is not over a monoid algebra of bimodules"))
    T = cell_tree(C)
    groups = _snapshot_groups(C, T)
    header = _snapshot_header(T, case, extra, groups, germ_fingerprint)
    stream = GzipCompressorStream(open(path, "w"))
    try
        for (k, v) in header
            write(stream, "# ", k, ": ", v, "\n")
        end
        for g in groups, l in g.atoms
            write(stream, l, "\n")
        end
    finally
        close(stream)
    end
    fingerprints && _write_fingerprints(string(chop(path; tail = length(".txt.gz")), ".sha.csv"), path, groups)
    return path
end

function _write_fingerprints(csv::AbstractString, gz::AbstractString, groups)
    rows = [(g.section, string(g.n), length(g.atoms), _sha_lines(g.atoms), _sha_lines(g.bar)) for g in groups]
    total_atoms = bytes2hex(SHA.sha256(join((r[4] for r in rows), "\n")))
    total_bar = bytes2hex(SHA.sha256(join((r[5] for r in rows), "\n")))
    file_sha = open(io -> bytes2hex(SHA.sha256(io)), gz)
    open(csv, "w") do io
        println(io, "section,n,lignes,sha256_atomes,sha256_barre")
        for r in rows
            println(io, join(r, ','))
        end
        println(io, "total,,", sum(r[3] for r in rows), ',', total_atoms, ',', total_bar)
        println(io, "fichier,,", filesize(gz), ',', file_sha, ',')
    end
    return csv
end
