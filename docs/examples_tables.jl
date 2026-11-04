# The tables of the pages Examples (docs/src/examples/), generated from the expected values
# of the tests, test/data/expected/*.toml, so that the pages and the tests share one source:
# a value added to these files appears on the pages at the next build. Included by
# docs/make.jl; the @eval blocks of the pages call these functions, which return Markdown.
# A copy in Python served to check the transformation of the data before delivery.
module ExamplesTables

import Markdown
import TOML

const DATA = joinpath(@__DIR__, "..", "test", "data", "expected")

# Subscript digits, for the names of the groups (G₃₄, A₁₂, I₂(5)) and the degrees (H₇).
subscript(s::AbstractString) = map(c -> isdigit(c) ? Char(0x2080 + (c - '0')) : c, s)

# The tables of a file, in the order of the file, which TOML.parsefile does not keep: the
# headers [G4], ["I2(5)"] and ["G(3,3,3)"], without [source] and the arrays [[provenance]].
function ordered_tables(name)
    path = joinpath(DATA, name * ".toml")
    data = TOML.parsefile(path)
    entries = Pair{String, Dict{String, Any}}[]
    for line in eachline(path)
        m = match(r"^\[\s*\"?([^\"\[\]]+)\"?\s*\]\s*$", line)
        m === nothing && continue
        key = String(m.captures[1])
        key != "source" && get(data, key, nothing) isa AbstractDict && push!(entries, key => data[key])
    end
    return entries
end

# "G4" → "G₄", "I2(5)" → "I₂(5)", "G(3,3,3)" → "G(3, 3, 3)".
function display_name(key)
    m = match(r"^G\((\d+),(\d+),(\d+)\)$", key)
    m === nothing || return "G(" * join(m.captures, ", ") * ")"
    m = match(r"^([A-Z])(\d+)\((\d+)\)$", key)
    m === nothing || return m[1] * subscript(m[2]) * "(" * m[3] * ")"
    m = match(r"^([A-Z])(\d+)$", key)
    m === nothing || return m[1] * subscript(m[2])
    return key
end

homology_label(n) = "H" * subscript(string(n))

# A Markdown table, the rows completed by empty cells.
function table(header, rows)
    w = length(header)
    lines = ["| " * join(header, " | ") * " |", "|" * repeat("---|", w)]
    for r in rows
        push!(lines, "| " * join(vcat(r, fill("", w - length(r))), " | ") * " |")
    end
    return join(lines, "\n")
end

maxlength(entries, key) = maximum(length(e[key]) for (_, e) in entries)

function homology_table(entries, key; columns = ["W"], extra = (k, e) -> String[], mark = (k, e) -> "")
    d = maxlength(entries, key)
    header = vcat(columns, [homology_label(n) for n in 0:(d - 1)])
    rows = [vcat([display_name(k) * mark(k, e)], extra(k, e), String.(e[key])) for (k, e) in entries]
    return table(header, rows)
end

function cells_table(entries)
    d = maxlength(entries, "cells_minimal")
    header = vcat(["W"], string.(0:(d - 1)), ["total"])
    rows = [vcat([display_name(k)], string.(e["cells_minimal"]), fill("", d - length(e["cells_minimal"])),
                 [string(sum(e["cells_minimal"]))]) for (k, e) in entries]
    return table(header, rows)
end

# A collapsed block of Documenter, around a Markdown text.
details(title, body) =
    "!!! details \"" * title * "\"\n" * join([isempty(l) ? "" : "    " * l for l in split(body, '\n')], "\n")

# Consecutive integers as ranges: [3, 4, 5, 6, 30] → "3–6, 30".
function ranges(ms)
    v = sort(collect(ms))
    runs = String[]
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[j + 1] == v[j] + 1
            j += 1
        end
        push!(runs, i == j ? string(v[i]) : j > i + 1 ? "$(v[i])–$(v[j])" : "$(v[i]), $(v[j])")
        i = j + 1
    end
    return join(runs, ", ")
end

# "A, B and C".
function enumeration(names)
    length(names) <= 1 && return join(names)
    return join(names[1:(end - 1)], ", ") * " and " * names[end]
end

# ------------------------------------------------------------------------------------
# Exceptional groups: the 26 well-generated ones, in the order of the tables of Garnier.

exceptional() = [k => e for (k, e) in ordered_tables("exceptional_braid_groups") if get(e, "well_generated", false)]

# The real type of W, or the representative of its isodiscriminantality class.
function exceptional_class(k, e)
    c = get(e, "real", get(e, "isodiscriminantal", ""))
    return [!isempty(c) && c != k ? display_name(c) : ""]
end

# The values with ℤε that the package does not compute for now (marked †).
const SIGN_NOT_COMPUTED = ("G37",)

exceptional_homology(key) =
    Markdown.parse(homology_table(exceptional(), key; columns = ["W", "Class"], extra = exceptional_class,
                                  mark = (k, e) -> key == "sign" && k in SIGN_NOT_COMPUTED ? " †" : ""))

exceptional_cells() = Markdown.parse(cells_table(exceptional()))

# Garnier, table 1: the Dehornoy–Lafont complexes of five groups beside the reduced resolution.
function exceptional_dehornoy_lafont()
    entries = [k => e for (k, e) in exceptional() if haskey(e, "cells_mar17")]
    d = maxlength(entries, "cells_minimal")
    labels = ("cells_minimal" => "reduced lcm resolution",
              "cells_mar17" => "Dehornoy–Lafont, order of Marin",
              "cells_optimised" => "Dehornoy–Lafont, order of Garnier")
    rows = [vcat([display_name(k), label], string.(e[key]), fill("", d - length(e[key])), [string(sum(e[key]))])
            for (k, e) in entries for (key, label) in labels]
    return Markdown.parse(table(vcat(["W", "complex"], string.(0:(d - 1)), ["total"]), rows))
end

# ------------------------------------------------------------------------------------
# Spherical Artin groups: the types A to H, then the dihedral types I₂(p).

spherical(dihedral) = [k => e for (k, e) in ordered_tables("spherical_artin_groups") if startswith(k, "I2(") == dihedral]

spherical_homology(key) = Markdown.parse(homology_table(spherical(false), key))

spherical_cells() = Markdown.parse(cells_table([k => e for (k, e) in spherical(false) if haskey(e, "cells_minimal")]))

function spherical_without_cells()
    names = [display_name(k) for (k, e) in spherical(false) if !haskey(e, "cells_minimal")]
    return Markdown.parse("The dual braid monoid of " * enumeration(names) *
                          " was not computed: their numbers of cells are not given.")
end

function spherical_dihedral()
    E = spherical(true)
    ps = [parse(Int, match(r"^I2\((\d+)\)$", k)[1]) for (k, _) in E]
    body = join(["Trivial coefficients:", "", homology_table(E, "trivial"), "",
                 "Sign module ℤε:", "", homology_table(E, "sign")], "\n")
    return Markdown.parse(details("I₂(p), p = " * ranges(ps), body))
end

# ------------------------------------------------------------------------------------
# Imprimitive groups: the G(m, m, n), one block for each n, and the G(m, 1, n).

mpn(k) = Tuple(parse.(Int, match(r"^G\((\d+),(\d+),(\d+)\)$", k).captures))

imprimitive() = [k => e for (k, e) in ordered_tables("imprimitive_braid_groups")]

imprimitive(n) = [k => e for (k, e) in imprimitive() if haskey(e, "trivial") && mpn(k)[3] == n && mpn(k)[1] == mpn(k)[2]]

imprimitive_ranks() = sort!(unique(mpn(k)[3] for (k, e) in imprimitive() if haskey(e, "trivial")))

unknown(e) = all(==("?"), e["sign"])

function imprimitive_block(n)
    E = imprimitive(n)
    known = [k => e for (k, e) in E if !unknown(e)]
    missing_m = [mpn(k)[1] for (k, e) in E if unknown(e)]
    parts = ["Trivial coefficients:", "", homology_table(E, "trivial"), "", "Sign module ℤε:", ""]
    isempty(known) || push!(parts, homology_table(known, "sign"), "")
    isempty(missing_m) || push!(parts, "Not computed yet for m = " * ranges(missing_m) * ".", "")
    push!(parts, "Numbers of cells:", "", cells_table(E))
    return details("G(m, m, $n), m = " * ranges(mpn(k)[1] for (k, _) in E), join(parts, "\n"))
end

imprimitive_blocks() = Markdown.parse(join([imprimitive_block(n) for n in imprimitive_ranks()], "\n\n"))

# The G(m, 1, n) whose dual germ is isomorphic to that of Bₙ.
function imprimitive_m1()
    E = [mpn(k) => e["dual_germ_isomorphic_to"] for (k, e) in imprimitive() if haskey(e, "dual_germ_isomorphic_to")]
    all(t == "B$(p[3])" for (p, t) in E) || error("a G(m, 1, n) not isomorphic to Bₙ in the data")
    ms = unique(p[1] for (p, _) in E)
    ns = unique(p[3] for (p, _) in E)
    return Markdown.parse("For m = " * ranges(ms) * " and n = " * ranges(ns) * ", the dual germ of " *
                          "G(m, 1, n) is isomorphic to that of the Coxeter group Bₙ (found by " *
                          "`GarsideHomology.germ_isomorphism`): the braid groups of G(m, 1, n) and of " *
                          "Bₙ are isomorphic, by an isomorphism that sends braid reflections to braid " *
                          "reflections, and their homology is given in the tables of the spherical " *
                          "Artin groups.")
end

# The values not computed yet with ℤε, for the section Not computed yet.
function imprimitive_unknown()
    lines = String[]
    for n in imprimitive_ranks()
        ms = [mpn(k)[1] for (k, e) in imprimitive(n) if unknown(e)]
        isempty(ms) || push!(lines, "- G(m, m, $n) for m = " * ranges(ms))
    end
    return Markdown.parse(join(lines, "\n"))
end

end
