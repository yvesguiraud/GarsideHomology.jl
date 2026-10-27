# Tests of the expected values of data/expected/ (f5 of J1.f): the three files checked
# without computation (format, alternating sums, links between the files, formulas), a
# small selection recomputed from the germs (2 or 3 cases per file), and, when the
# environment variable GAUSS_LONG_TESTS is set, every entry recomputed.
#
# GAUSS_LONG_TESTS = "all" recomputes every entry, any other value is a regular expression
# on the names of the tables ('^(B[2-4]|G\(3,3,.\))$', say). The long tests need a solver
# that widens its integers (SparseElementaryDivisors 0.1.1, where DEFAULT_VALUE_TYPES is
# defined); with an older one they are skipped. Until the version 0.1.1 is registered, they
# run in an environment where the solver is developed, such as CatLab/j1:
#     julia --project=j1 -t 16 <path of Gauss.jl>/test/test_expected_values.jl
# About 1 h 45 min for all the entries on a machine like N9 (estimate from J1.f).

using Gauss
using Gauss: rank, germ_isomorphism, cell_tree, ncells
using Test
import Chevie
import SparseElementaryDivisors
import TOML

const EV_DIR = joinpath(@__DIR__, "data", "expected")
const EV_G = FinitelyGeneratedAbelianGroup
const EV_GG = FinitelyGeneratedGradedAbelianGroup
const EV_MODULES = ("trivial" => TrivialModule(), "sign" => SignModule())

ev_read(file) = TOML.parsefile(joinpath(EV_DIR, file))
ev_tables(d) = sort!([k for k in keys(d) if !(k in ("source", "provenance"))])
ev_graded(v) = EV_GG([parse(EV_G, s) for s in v])
ev_displayed(v) = all(s == "?" || string(parse(EV_G, s)) == s for s in v)
ev_euler(v) = sum((-1)^k * rank(parse(EV_G, v[k + 1])) for k in 0:(length(v) - 1)) == 0
ev_euler(v::Vector{<:Integer}) = sum((-1)^k * v[k + 1] for k in 0:(length(v) - 1)) == 0

# The rank of a Coxeter type named as in spherical_artin_groups.toml ("B3", "I2(5)") and
# its Coxeter group.
function ev_coxeter(name)
    m = match(r"^I2\((\d+)\)$", name)
    m === nothing || return 2, Chevie.coxgroup(:I, 2, parse(Int, m[1]))
    m = match(r"^([A-H])(\d+)$", name)
    n = parse(Int, m[2])
    return n, Chevie.coxgroup(Symbol(m[1]), n)
end

# (m, p, n) of a table "G(m,p,n)".
ev_mpn(name) = Tuple(parse.(Int, match(r"^G\((\d+),(\d+),(\d+)\)$", name).captures))

# Callegaro–Marin, theorem 6.4: H₂ of B(e, e, r), r ≥ 3, by its invariant factors; Betti
# numbers of the Poincaré polynomial 1 + t, or 1 + t + tʳ⁻¹ + tʳ when e and r are even.
function ev_callegaro_marin_h2(e, r)
    r == 3 && return EV_G(0, [e])
    isodd(e) && return EV_G(0, [2e])
    return EV_G(0, r == 4 ? [2, 2, e] : [2, e])
end
ev_betti(e, r, k) = k in (0, 1) ? 1 : (iseven(e) && iseven(r) && k in (r - 1, r)) ? 1 : 0

# The cells by degree and the homologies of a germ, against an entry. `skip_sign`: ℤε not
# recomputed (a value that needs more than the integers allowed, or not computed in J1).
function ev_check_germ(G, entry; value_type, cells = get(entry, "cells_minimal", nothing),
                       skip_sign = false)
    R = lcm_resolution(G)
    T = cell_tree(R)
    cells === nothing || @test [ncells(T, n) for n in range(R)] == cells
    for (key, M) in EV_MODULES
        ((key == "sign" && skip_sign) || "?" in entry[key]) && continue
        @test homology(R, M; value_type) == ev_graded(entry[key])
    end
    return nothing
end

const EV_EXCEPTIONAL = ev_read("exceptional_braid_groups.toml")
const EV_SPHERICAL = ev_read("spherical_artin_groups.toml")
const EV_IMPRIMITIVE = ev_read("imprimitive_braid_groups.toml")
const EV_PROVENANCE_KEYS = Set(["run", "description", "script", "data", "catlab", "date", "machine", "threads",
                                "julia", "gauss", "solver", "chevie", "numbering", "value_types", "ceiling"])

@testset "Expected values" begin

    @testset "Provenance" begin
        for (d, n) in ((EV_EXCEPTIONAL, 1), (EV_SPHERICAL, 3), (EV_IMPRIMITIVE, 2))
            runs = d["provenance"]
            @test runs isa Vector && length(runs) == n
            for run in runs
                @test issubset(keys(run), EV_PROVENANCE_KEYS)
                @test all(haskey(run, k) for k in ("run", "description", "catlab", "date", "machine", "gauss", "julia"))
            end
        end
    end

    @testset "Exceptional groups" begin
        # Tables 1, 3 and 4 of the article named in [source], with the errata (f1 of J1.f),
        # and cells_minimal (f5).
        expected = EV_EXCEPTIONAL
        @test expected["source"]["arxiv"] == "2302.10033v2"
        @test expected["source"]["tables"] == [1, 3, 4]
        groups = ev_tables(expected)
        @test sort!(parse.(Int, chop.(groups; head = 1, tail = 0))) == 4:37
        @test count(g -> expected[g]["well_generated"], groups) == 26
        keys_known = Set(["well_generated", "real", "isodiscriminantal", "trivial", "sign", "cells_minimal",
                          "cells_mar17", "cells_optimised", "two_cell_bounds", "errata"])
        classes = Dict{String, Vector{String}}()
        n_errata = 0
        for group in groups
            entry = expected[group]
            @test issubset(keys(entry), keys_known)
            # Written in the form displayed by FinitelyGeneratedAbelianGroup.
            @test ev_displayed(entry["trivial"]) && ev_displayed(entry["sign"])
            trivial = ev_graded(entry["trivial"])
            sign = ev_graded(entry["sign"])
            @test range(sign) == range(trivial)
            @test trivial[0] == EV_G(1)
            @test sign[0] == EV_G(0, [2])
            # Braid groups of irreducible reflection groups: Euler characteristic 0, the
            # alternating sum of the free ranks for any module of rank one.
            @test ev_euler(entry["trivial"]) && ev_euler(entry["sign"])
            for key in ("cells_mar17", "cells_optimised", "cells_minimal")
                haskey(entry, key) || continue
                cells = entry[key]
                @test length(cells) == length(range(trivial))
                @test first(cells) == 1 && ev_euler(cells)
            end
            # cells_minimal for the well-generated groups exactly, never more than the cells
            # of a Dehornoy–Lafont complex (it is minimal).
            @test haskey(entry, "cells_minimal") == entry["well_generated"]
            haskey(entry, "cells_optimised") && @test all(entry["cells_minimal"] .<= entry["cells_optimised"])
            for erratum in get(entry, "errata", Dict{String, Any}[])
                @test erratum["key"] in ("trivial", "sign", "cells_mar17", "cells_optimised")
                @test string(entry[erratum["key"]][erratum["degree"] + 1]) != erratum["printed"]
                @test erratum["listed"] isa Bool && !isempty(erratum["reason"])
                n_errata += 1
            end
            haskey(entry, "isodiscriminantal") &&
                push!(get!(classes, entry["isodiscriminantal"], String[]), group)
        end
        @test n_errata == 3
        @test length(classes) == 10
        # The groups of an isodiscriminantality class have the same homology, that of the
        # representative, a table of this file or of spherical_artin_groups.toml.
        for (representative, members) in classes, group in members
            reference = get(expected, representative, get(EV_SPHERICAL, representative, nothing))
            @test reference isa Dict
            @test expected[group]["trivial"] == reference["trivial"]
            @test expected[group]["sign"] == reference["sign"]
        end
        # A real group is its Coxeter type: same homology, same dual germ up to isomorphism.
        for group in groups
            haskey(expected[group], "real") || continue
            reference = EV_SPHERICAL[expected[group]["real"]]
            for key in ("trivial", "sign", "cells_minimal")
                @test expected[group][key] == reference[key]
            end
        end
        @test expected["G37"]["real"] == "E8"
    end

    @testset "Artin groups of spherical type" begin
        names = vcat(["A$n" for n in 1:12], ["B$n" for n in 2:10], ["D$n" for n in 4:10],
                     ["E6", "E7", "E8", "F4", "H3", "H4"], ["I2($p)" for p in 5:30])
        @test ev_tables(EV_SPHERICAL) == sort(names)
        @test sort!([n for n in names if !haskey(EV_SPHERICAL[n], "cells_minimal")]) ==
              sort!(["A10", "A11", "A12", "B9", "B10", "D9", "D10"])
        for name in names
            entry = EV_SPHERICAL[name]
            n, _ = ev_coxeter(name)
            @test issubset(keys(entry), Set(["trivial", "sign", "cells_minimal"]))
            @test length(entry["trivial"]) == length(entry["sign"]) == n + 1
            @test ev_displayed(entry["trivial"]) && ev_displayed(entry["sign"])
            @test !("?" in entry["trivial"]) && !("?" in entry["sign"])
            @test entry["trivial"][1] == "ℤ" && entry["sign"][1] == "ℤ/2ℤ"
            @test ev_euler(entry["trivial"]) && ev_euler(entry["sign"])
            if haskey(entry, "cells_minimal")
                cells = entry["cells_minimal"]
                @test length(cells) == n + 1 && first(cells) == 1 && ev_euler(cells)
            end
        end
        # Rank 2: in ℤ, H₁ = ℤ and H₂ = 0 for p odd, ℤ² and ℤ for p even; in ℤε, ℤ/2ℤ, ℤ/pℤ,
        # 0 (observed in J1.f for p = 3 to 30); p atoms and p − 1 cells of degree 2.
        for p in 5:30
            entry = EV_SPHERICAL["I2($p)"]
            @test entry["trivial"] == (isodd(p) ? ["ℤ", "ℤ", "0"] : ["ℤ", "ℤ²", "ℤ"])
            @test entry["sign"] == ["ℤ/2ℤ", "ℤ/$(p)ℤ", "0"]
            @test entry["cells_minimal"] == [1, p, p - 1]
        end
        @test EV_SPHERICAL["A2"]["sign"] == ["ℤ/2ℤ", "ℤ/3ℤ", "0"] && EV_SPHERICAL["B2"]["sign"] == ["ℤ/2ℤ", "ℤ/4ℤ", "0"]
    end

    @testset "Imprimitive groups" begin
        names = ev_tables(EV_IMPRIMITIVE)
        m1n = filter(name -> ev_mpn(name)[2] == 1, names)
        mmn = filter(name -> ev_mpn(name)[1] == ev_mpn(name)[2], names)
        @test length(m1n) == 42 && length(mmn) == 144 && length(names) == 186
        @test sort!(m1n) == sort!(["G($m,1,$n)" for n in 2:8 for m in 3:8])
        for name in m1n
            m, _, n = ev_mpn(name)
            @test EV_IMPRIMITIVE[name] == Dict("dual_germ_isomorphic_to" => "B$n")
            @test haskey(EV_SPHERICAL, "B$n")
        end
        n_unknown = 0
        for name in mmn
            entry = EV_IMPRIMITIVE[name]
            e, _, r = ev_mpn(name)
            @test 3 <= r <= 8 && 3 <= e <= 30
            @test keys(entry) == Set(["trivial", "sign", "cells_minimal"])
            trivial, sign, cells = entry["trivial"], entry["sign"], entry["cells_minimal"]
            @test length(trivial) == length(sign) == length(cells) == r + 1
            @test ev_displayed(trivial) && ev_displayed(sign)
            # In ℤ: Callegaro–Marin, theorem 6.4, and the Poincaré polynomial.
            @test !("?" in trivial) && ev_euler(trivial)
            @test trivial[2] == "ℤ"
            @test parse(EV_G, trivial[3]) == ev_callegaro_marin_h2(e, r)
            @test all(rank(parse(EV_G, trivial[k + 1])) == ev_betti(e, r, k) for k in 0:r)
            # In ℤε: computed, or "?" in every degree (beyond Int128 in J1.f).
            if "?" in sign
                @test all(==("?"), sign)
                n_unknown += 1
            else
                @test sign[1] == "ℤ/2ℤ" && ev_euler(sign)
            end
            # Cells: the atoms (n − 1)(m + n − 2) and |μ(1, Δ)| = C(2n−3, n−1)·(2(n−1)m − n)/n.
            μ, remainder = divrem(binomial(2r - 3, r - 1) * (2 * (r - 1) * e - r), r)
            @test remainder == 0
            @test cells[1] == 1 && cells[2] == (r - 1) * (e + r - 2) && cells[end] == μ && ev_euler(cells)
        end
        @test n_unknown == 29
    end

    # A small selection recomputed: rank 2, 3 and 4 exceptional groups (G32: the erratum on
    # H₄), the classical and dual monoids of B₃ and D₄ (no published value), I₂(5), G(e, e, r)
    # for e and r odd, then even (the longer Poincaré polynomial), and G(3, 1, 3) against
    # B₃. ℤε with value_type = Int64, accepted by SparseElementaryDivisors 0.1.0 and 0.1.1
    # (the version 0.1.0 starts in Int32 and does not widen; all these cases hold in Int64).
    # The isomorphisms of germs between Chevie's presentations (crg(1, 1, 3) and A₂, crg(2,
    # 2, 4) and D₄, ...) are tested in test_chevie_germs.jl.
    @testset "Recomputed: $name" for name in ("G4", "G27", "G32")
        ev_check_germ(dual_germ(Chevie.crg(parse(Int, name[2:end]))), EV_EXCEPTIONAL[name]; value_type = Int64)
    end
    @testset "Recomputed: $name" for name in ("B3", "D4", "I2(5)")
        n, W = ev_coxeter(name)
        entry = EV_SPHERICAL[name]
        ev_check_germ(dual_germ(W), entry; value_type = Int64)
        name == "I2(5)" ||
            ev_check_germ(germ(Chevie.BraidMonoid(W)), entry; value_type = Int64, cells = [binomial(n, k) for k in 0:n])
    end
    @testset "Recomputed: $name" for name in ("G(3,3,3)", "G(4,4,4)")
        ev_check_germ(dual_germ(Chevie.crg(ev_mpn(name)...)), EV_IMPRIMITIVE[name]; value_type = Int64)
    end
    @testset "Recomputed: G(3,1,3)" begin
        type = EV_IMPRIMITIVE["G(3,1,3)"]["dual_germ_isomorphic_to"]
        @test germ_isomorphism(dual_germ(Chevie.crg(3, 1, 3)), dual_germ(last(ev_coxeter(type)))) !== nothing
    end
end

# The long tests: every entry recomputed (or those whose name matches GAUSS_LONG_TESTS),
# integers widened up to Int128, as in J1.f. Not recomputed: ℤε where it was not computed
# in J1 ("?", B*(A₉), B*(E₈), whose `sign` comes from B⁺(W), and G37, from Magma); the
# groups that are not well generated (no dual braid monoid). For G(m, 1, n), the
# isomorphism of dual germs with Bₙ, which gives the same resolution up to the numbering.
const EV_LONG = get(ENV, "GAUSS_LONG_TESTS", "")
const EV_LONG_TYPES = (Int32, Int64, Int128)
const EV_SIGN_NOT_RECOMPUTED = ("G37", "dual A9", "dual E8")

function ev_long_line(label, t)
    println("  ", rpad(label, 22), lpad(round(t; digits = 1), 8), " s")
    flush(stdout)
end

if isempty(EV_LONG)
    @info "Long tests of the expected values skipped: set GAUSS_LONG_TESTS to all, or to a " *
          "regular expression on the names of the tables."
elseif !isdefined(SparseElementaryDivisors, :DEFAULT_VALUE_TYPES)
    @info "Long tests of the expected values skipped: the solver does not widen its integers " *
          "(SparseElementaryDivisors 0.1.1 needed)."
else
    @testset "Expected values, long tests" begin
        pattern = EV_LONG == "all" ? r"" : Regex(EV_LONG)
        selected(d) = [name for name in ev_tables(d) if occursin(pattern, name)]
        println("Long tests of the expected values, for the tables matching $(repr(EV_LONG)):")
        elapsed = @elapsed begin
            @testset "Long: $name" for name in filter(g -> EV_EXCEPTIONAL[g]["well_generated"], selected(EV_EXCEPTIONAL))
                entry = EV_EXCEPTIONAL[name]
                t = @elapsed ev_check_germ(dual_germ(Chevie.crg(parse(Int, name[2:end]))), entry;
                                           value_type = EV_LONG_TYPES, skip_sign = name in EV_SIGN_NOT_RECOMPUTED)
                ev_long_line(name, t)
            end
            @testset "Long: $name" for name in selected(EV_SPHERICAL)
                n, W = ev_coxeter(name)
                entry = EV_SPHERICAL[name]
                if !startswith(name, "I2")
                    t = @elapsed ev_check_germ(germ(Chevie.BraidMonoid(W)), entry; value_type = EV_LONG_TYPES,
                                               cells = [binomial(n, k) for k in 0:n])
                    ev_long_line("classic $name", t)
                end
                if haskey(entry, "cells_minimal")
                    t = @elapsed ev_check_germ(dual_germ(W), entry; value_type = EV_LONG_TYPES,
                                               skip_sign = "dual $name" in EV_SIGN_NOT_RECOMPUTED)
                    ev_long_line("dual $name", t)
                end
            end
            @testset "Long: $name" for name in selected(EV_IMPRIMITIVE)
                entry = EV_IMPRIMITIVE[name]
                W = Chevie.crg(ev_mpn(name)...)
                if haskey(entry, "dual_germ_isomorphic_to")
                    t = @elapsed @test germ_isomorphism(dual_germ(W), dual_germ(last(ev_coxeter(entry["dual_germ_isomorphic_to"])))) !== nothing
                else
                    t = @elapsed ev_check_germ(dual_germ(W), entry; value_type = EV_LONG_TYPES)
                end
                ev_long_line(name, t)
            end
        end
        println("Long tests: $(round(elapsed / 60; digits = 1)) min.")
    end
end
