# Reading and writing the reference snapshots of J0 (format version 1, see
# data/reference/FORMAT.md) with the types of Gauss.Complexes.
#
# A tool for the tests and for the comparisons of the new construction (J1.d), not part
# of the package. A snapshot is read into a `Snapshot`: its header, the simples met (by
# their words) and a chain complex of based free modules over the monoid algebra, whose
# differentials are labelled by pairs (σ, τ) of simples. Writing regenerates every data
# line from these types, so that equal fingerprints show that nothing was lost.

using Gauss
using Gauss: MonoidAlgebra    # public, not exported (Oscar exports a MonoidAlgebra)
using CodecZlib: GzipCompressorStream, GzipDecompressorStream
import SHA

const _SNAP_ZZ = Gauss.Complexes.AbstractAlgebra.ZZ
const _SNAP_SA = Gauss.Complexes.SparseArrays

"""
    SnapshotMonoid(simples)

The simples met in a snapshot, given by their smallest words in the atoms (`"2.1.3"`,
`"e"` for the identity) and numbered 0 for the identity, then by length and word:
the simple numbered `i` is `simples[i + 1]`.
"""
struct SnapshotMonoid <: Gauss.Complexes.AbstractMonoid
    simples::Vector{String}
    index::Dict{String, Int32}
end

SnapshotMonoid(simples::Vector{String}) =
    SnapshotMonoid(simples, Dict(w => Int32(i - 1) for (i, w) in enumerate(simples)))

Base.:(==)(a::SnapshotMonoid, b::SnapshotMonoid) = a.simples == b.simples
Base.hash(m::SnapshotMonoid, h::UInt) = hash(m.simples, hash(:SnapshotMonoid, h))

"Length of the simple numbered `i` of a snapshot: the number of atoms of its word."
simple_length(m::SnapshotMonoid, i::Integer) =
    m.simples[i + 1] == "e" ? 0 : count(==('.'), m.simples[i + 1]) + 1

# Order of the simples: by length, then by word, read as a sequence of atom indices.
_word_key(w::AbstractString) =
    w == "e" ? (0, Int[]) : (count(==('.'), w) + 1, parse.(Int, split(w, '.')))

"""
    SnapshotCell(atoms, bar, weight)

A cell of a snapshot: its atom sequence (`"7.3.12"`, `"-"` in degree 0), its bar form
(`"[1.2|3|2.1]"`, `"[]"` in degree 0) and the number of its weight in the monoid.
"""
struct SnapshotCell
    atoms::String
    bar::String
    weight::Int32
end

const SnapshotAlgebra = MonoidAlgebra{typeof(_SNAP_ZZ), SnapshotMonoid}
const SnapshotModule = BasedFreeModule{SnapshotAlgebra, Vector{SnapshotCell}}
const SnapshotMatrix = LabelledSparseMatrix{NTuple{2, Int32}, Int}
const SnapshotComplex = ChainComplex{SnapshotAlgebra, SnapshotModule, SnapshotMatrix}

"""
    Snapshot

A snapshot read from a file: `header` (the pairs `key => value` of the header, in their
order), `monoid` (a [`SnapshotMonoid`](@ref)) and `complex`, the chain complex over
`MonoidAlgebra(ZZ, monoid, :both)` whose module of degree n has the cells of degree n as
basis, in the order of the file, and whose `differential(complex, n)` holds the terms of
the lines `D n`, with labels `(σ, τ)`.
"""
struct Snapshot
    header::Vector{Pair{String, String}}
    monoid::SnapshotMonoid
    complex::SnapshotComplex
end

"The value of `key` in the header of `s`, or `nothing`."
function header_value(s::Snapshot, key::AbstractString)
    k = findfirst(p -> p.first == key, s.header)
    return k === nothing ? nothing : s.header[k].second
end

"""
    snapshots_equal(a::Snapshot, b::Snapshot)

Whether two snapshots have the same header, simples, cells and differentials.
"""
function snapshots_equal(a::Snapshot, b::Snapshot)
    a.header == b.header && a.monoid == b.monoid || return false
    A, B = a.complex, b.complex
    range(A) == range(B) || return false
    all(cells(A[n]) == cells(B[n]) for n in range(A)) || return false
    return all(differential(A, n) == differential(B, n) for n in first(range(A)):(last(range(A)) - 1))
end

# ------------------------------------------------------------------------------------
# Reading

_snapshot_error(name, lineno, msg) =
    throw(ArgumentError("snapshot $name, line $lineno: $msg"))

# Number of the simple `w` in the order met, numbering it if it is new.
function _simple_number!(met::Dict{String, Int32}, words::Vector{String}, w::AbstractString)
    k = get(met, w, Int32(0))
    k == 0 || return k
    push!(words, String(w))
    k = Int32(length(words))
    met[words[end]] = k
    return k
end

"""
    read_snapshot(path) -> Snapshot

Read a gzip-compressed snapshot of format version 1. Checks the format line, the order
of the groups, the degrees of the source and target of each term, that cells and terms
are not repeated, that coefficients are not zero, and the numbers of cells and terms
announced by the header; any failure throws an `ArgumentError` with the line number.
"""
function read_snapshot(path::AbstractString)
    stream = GzipDecompressorStream(open(path))
    try
        return read_snapshot(stream, basename(path))
    finally
        close(stream)
    end
end

function read_snapshot(io::IO, name::AbstractString)
    header = Pair{String, String}[]
    words = String[]                            # simples, in the order met
    met = Dict{String, Int32}()                 # word -> provisional number
    cells = Vector{SnapshotCell}[]              # cells[n + 1]: the cells of degree n
    position = Dict{String, Int32}[]            # atom sequence -> position, per degree
    I = Vector{Int32}[]                         # terms of d_n, in I[n + 1], J[n + 1], ...
    J = Vector{Int32}[]
    L = Vector{NTuple{2, Int32}}[]
    V = Vector{Int}[]
    phase = 0                                   # 0 header, 1 cells, 2 terms
    last_n = -1
    nterms = 0
    lineno = 0
    for line in eachline(io)
        lineno += 1
        if startswith(line, "# ")
            phase == 0 || _snapshot_error(name, lineno, "header line after the data")
            kv = split(SubString(line, 3), ": "; limit = 2)
            length(kv) == 2 || _snapshot_error(name, lineno, "header line without \": \"")
            push!(header, String(kv[1]) => String(kv[2]))
            continue
        end
        f = split(line, ' ')
        if f[1] == "C"
            phase <= 1 || _snapshot_error(name, lineno, "cell after the terms")
            phase == 0 && (isempty(header) || header[1] != ("format" => "cliche-gauss 1")) &&
                _snapshot_error(name, lineno, "the first header line must be \"format: cliche-gauss 1\"")
            phase = 1
            length(f) == 5 || _snapshot_error(name, lineno, "a cell line has 5 fields")
            n = parse(Int, f[2])
            if n == length(cells)
                push!(cells, SnapshotCell[])
                push!(position, Dict{String, Int32}())
            elseif n != length(cells) - 1
                _snapshot_error(name, lineno, "cell of degree $n after degree $(length(cells) - 1)")
            end
            haskey(position[n + 1], f[3]) && _snapshot_error(name, lineno, "repeated cell $(f[3])")
            w = _simple_number!(met, words, f[5])
            push!(cells[n + 1], SnapshotCell(String(f[3]), String(f[4]), w))
            position[n + 1][String(f[3])] = Int32(length(cells[n + 1]))
        elseif f[1] == "D"
            phase >= 1 || _snapshot_error(name, lineno, "term before the cells")
            if phase == 1
                phase = 2
                for _ in 1:(length(cells) - 1)
                    push!(I, Int32[]); push!(J, Int32[]); push!(L, NTuple{2, Int32}[]); push!(V, Int[])
                end
            end
            length(f) == 7 || _snapshot_error(name, lineno, "a term line has 7 fields")
            n = parse(Int, f[2])
            0 <= n < length(cells) - 1 ||
                _snapshot_error(name, lineno, "term of d_$n, degrees of cells 0 to $(length(cells) - 1)")
            n >= last_n || _snapshot_error(name, lineno, "term of d_$n after d_$last_n")
            last_n = n
            j = get(position[n + 2], f[3], Int32(0))
            j == 0 && _snapshot_error(name, lineno, "unknown source $(f[3]) in degree $(n + 1)")
            i = get(position[n + 1], f[6], Int32(0))
            i == 0 && _snapshot_error(name, lineno, "unknown target $(f[6]) in degree $n")
            c = parse(Int, f[4])
            c == 0 && _snapshot_error(name, lineno, "zero coefficient")
            σ = _simple_number!(met, words, f[5])
            τ = _simple_number!(met, words, f[7])
            push!(I[n + 1], i); push!(J[n + 1], j); push!(L[n + 1], (σ, τ)); push!(V[n + 1], c)
            nterms += 1
        else
            _snapshot_error(name, lineno, "unknown line type \"$(f[1])\"")
        end
    end
    isempty(cells) && _snapshot_error(name, lineno, "no cell")
    for _ in (length(I) + 1):(length(cells) - 1)        # no term at all
        push!(I, Int32[]); push!(J, Int32[]); push!(L, NTuple{2, Int32}[]); push!(V, Int[])
    end
    ncells = sum(length, cells)
    for (key, value) in (("cellules", ncells), ("termes", nterms))
        k = findfirst(p -> p.first == key, header)
        k === nothing && _snapshot_error(name, lineno, "no \"$key\" in the header")
        parse(Int, header[k].second) == value ||
            _snapshot_error(name, lineno, "$value $key read, $(header[k].second) announced")
    end

    # Final numbering of the simples: identity 0, then by length and word.
    haskey(met, "e") || _snapshot_error(name, lineno, "the identity e does not appear")
    order = sortperm(words; by = _word_key)
    renumber = Vector{Int32}(undef, length(words))
    for (r, p) in enumerate(order)
        renumber[p] = r - 1
    end
    monoid = SnapshotMonoid(words[order])
    A = MonoidAlgebra(_SNAP_ZZ, monoid, :both)
    modules = [BasedFreeModule(A, [SnapshotCell(x.atoms, x.bar, renumber[x.weight]) for x in cs])
               for cs in cells]
    maps = SnapshotMatrix[]
    for n in 0:(length(cells) - 2)
        labs = [(renumber[σ], renumber[τ]) for (σ, τ) in L[n + 1]]
        d = labelled_sparse(I[n + 1], J[n + 1], labs, V[n + 1], length(cells[n + 1]),
                            length(cells[n + 2]))
        _SNAP_SA.nnz(d) == length(V[n + 1]) ||
            _snapshot_error(name, lineno, "repeated terms in d_$n")
        push!(maps, d)
    end
    return Snapshot(header, monoid, ChainComplex(modules, maps))
end

# ------------------------------------------------------------------------------------
# Lines, fingerprints, writing

_coefficient_string(c) = c > 0 ? "+$c" : string(c)

"""
    foreach_group(f, s::Snapshot)

Call `f(section, n, atom_lines, bar_lines)` for each group of data lines of `s`, in the
order of the file (`"C"` for each degree, then `"D"` for each `d_n`), with the lines of
the two forms of FORMAT.md, each sorted by bytes. The lines are regenerated from the
cells and the labelled differentials.
"""
function foreach_group(f, s::Snapshot)
    C = s.complex
    w = s.monoid.simples
    for n in range(C)
        cs = cells(C[n])
        atom_lines = ["C $n $(x.atoms) $(x.bar) $(w[x.weight + 1])" for x in cs]
        bar_lines = ["C $n $(x.bar) $(w[x.weight + 1])" for x in cs]
        f("C", n, sort!(atom_lines), sort!(bar_lines))
    end
    for n in first(range(C)):(last(range(C)) - 1)
        d = differential(C, n)
        source, target = cells(C[n + 1]), cells(C[n])
        rows, labs, vals = _SNAP_SA.rowvals(d), labels(d), _SNAP_SA.nonzeros(d)
        atom_lines = Vector{String}(undef, _SNAP_SA.nnz(d))
        bar_lines = Vector{String}(undef, _SNAP_SA.nnz(d))
        for j in 1:ncols(d), k in _SNAP_SA.nzrange(d, j)
            x, y = source[j], target[rows[k]]
            σ, τ = w[labs[k][1] + 1], w[labs[k][2] + 1]
            c = _coefficient_string(vals[k])
            atom_lines[k] = "D $n $(x.atoms) $c $σ $(y.atoms) $τ"
            bar_lines[k] = "D $n $(x.bar) $c $σ $(y.bar) $τ"
        end
        f("D", n, sort!(atom_lines), sort!(bar_lines))
    end
    return nothing
end

const _NEWLINE = UInt8[0x0a]

function _sha_lines(lines::Vector{String})
    ctx = SHA.SHA256_CTX()
    for l in lines
        SHA.update!(ctx, codeunits(l))
        SHA.update!(ctx, _NEWLINE)
    end
    return bytes2hex(SHA.digest!(ctx))
end

"""
    snapshot_fingerprints(s::Snapshot)

The rows of the `.sha.csv` file of `s` (FORMAT.md, section 5), without the row
`fichier`: named tuples `(section, n, lines, atoms, bar)`, `n` and the SHA-256 being
strings, ending with the row `total`.
"""
function snapshot_fingerprints(s::Snapshot)
    rows = NamedTuple{(:section, :n, :lines, :atoms, :bar), Tuple{String, String, Int, String, String}}[]
    foreach_group(s) do section, n, atom_lines, bar_lines
        push!(rows, (section = section, n = string(n), lines = length(atom_lines),
                     atoms = _sha_lines(atom_lines), bar = _sha_lines(bar_lines)))
    end
    push!(rows, (section = "total", n = "", lines = sum(r.lines for r in rows),
                 atoms = bytes2hex(SHA.sha256(join((r.atoms for r in rows), "\n"))),
                 bar = bytes2hex(SHA.sha256(join((r.bar for r in rows), "\n")))))
    return rows
end

"""
    read_fingerprints(path) -> (rows, file)

The rows of a `.sha.csv` file in the form of [`snapshot_fingerprints`](@ref), and its
row `fichier` as `file = (size, sha)`.
"""
function read_fingerprints(path::AbstractString)
    lines = readlines(path)
    lines[1] == "section,n,lignes,sha256_atomes,sha256_barre" ||
        throw(ArgumentError("$path: unexpected header $(lines[1])"))
    rows = NamedTuple{(:section, :n, :lines, :atoms, :bar), Tuple{String, String, Int, String, String}}[]
    file = (size = -1, sha = "")
    for l in lines[2:end]
        isempty(l) && continue
        f = split(l, ',')
        length(f) == 5 || throw(ArgumentError("$path: row with $(length(f)) fields: $l"))
        if f[1] == "fichier"
            file = (size = parse(Int, f[3]), sha = String(f[4]))
        else
            push!(rows, (section = String(f[1]), n = String(f[2]), lines = parse(Int, f[3]),
                         atoms = String(f[4]), bar = String(f[5])))
        end
    end
    return rows, file
end

"""
    write_snapshot(path, s::Snapshot)

Write `s` gzip-compressed in format version 1: its header as read, then the data lines
regenerated by [`foreach_group`](@ref).
"""
function write_snapshot(path::AbstractString, s::Snapshot)
    stream = GzipCompressorStream(open(path, "w"))
    try
        for (key, value) in s.header
            write(stream, "# ", key, ": ", value, "\n")
        end
        foreach_group(s) do section, n, atom_lines, bar_lines
            for l in atom_lines
                write(stream, l, "\n")
            end
        end
    finally
        close(stream)
    end
    return path
end

"The decompressed text of a gzip-compressed file."
gunzip_text(path::AbstractString) =
    open(io -> read(GzipDecompressorStream(io), String), path)

# ------------------------------------------------------------------------------------
# Collapsible pairs, as counted by the script B3 of J0

"""
    identity_monomials(C)

The monomials of the differentials of a complex of snapshot type whose label is the
identity (σ = τ = e), as the script B3 of J0 counted them, without looking at the
coefficient: named tuples `(n, source, target, coefficient)`.
"""
function identity_monomials(C)
    monomials = NamedTuple{(:n, :source, :target, :coefficient), NTuple{4, Int}}[]
    one_label = Gauss.Complexes.identity_label(NTuple{2, Int32})
    for n in first(range(C)):(last(range(C)) - 1)
        d = differential(C, n)
        rows, labs, vals = _SNAP_SA.rowvals(d), labels(d), _SNAP_SA.nonzeros(d)
        for j in 1:ncols(d), k in _SNAP_SA.nzrange(d, j)
            labs[k] == one_label &&
                push!(monomials, (n = n, source = j, target = Int(rows[k]), coefficient = vals[k]))
        end
    end
    return monomials
end

"""
    pair_counts(pairs)

The counts `(sources, pairs, targets)` of a list of pairs or monomials, as B3 gave them:
the distinct `(n, source)`, the elements of the list, the distinct `(n, target)`.
"""
pair_counts(pairs) = (length(unique((p.n, p.source) for p in pairs)), length(pairs),
                      length(unique((p.n, p.target) for p in pairs)))

"""
    read_pair_counts(path)

The file `collapsible_pairs.csv` of data/reference as a dictionary
`"<case>_<variant>" => (sources, pairs, targets)`, `targets` being `nothing` for the
snapshots recorded only by the first series of J0, whose output had no targets.
"""
function read_pair_counts(path::AbstractString)
    lines = readlines(path)
    lines[1] == "case,variant,sources,pairs,targets" ||
        throw(ArgumentError("$path: unexpected header $(lines[1])"))
    counts = Dict{String, NamedTuple{(:sources, :pairs, :targets), Tuple{Int, Int, Union{Int, Nothing}}}}()
    for l in lines[2:end]
        isempty(l) && continue
        f = split(l, ',')
        length(f) == 5 || throw(ArgumentError("$path: row with $(length(f)) fields: $l"))
        counts[f[1] * "_" * f[2]] = (sources = parse(Int, f[3]), pairs = parse(Int, f[4]),
                                     targets = isempty(f[5]) ? nothing : parse(Int, f[5]))
    end
    return counts
end

"""
    pair_counts_agree(counts, reference)

Whether the counts `(sources, pairs, targets)` agree with a row of
[`read_pair_counts`](@ref), the targets being compared only when they are known.
"""
pair_counts_agree(counts::NTuple{3, Int}, reference) =
    counts[1] == reference.sources && counts[2] == reference.pairs &&
    (reference.targets === nothing || counts[3] == reference.targets)
