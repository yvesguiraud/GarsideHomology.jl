# Tests of the projection onto the cells (cells/chains.jl, cells/projection.jl) and of the
# D lines of the snapshots: the tree of the chains, the splitting, the projection and the
# boundary on B⁺(A₂) and B*(A₂) by hand, ∂² = 0, and the fingerprints of the D lines of
# small cases of J0 (strict for Gauss_a and Gauss_d, bar form for Gauss_c, ∂² = 0 only for
# Gauss_b and for Gauss_a of the classical monoids).

isdefined(Main, :TREE_VARIANTS) || include("test_tree.jl")

using Gauss
using Gauss: CellTree, ncells, ChainTree, BimoduleTerm, LCMProjectionData,
             lcm_splitting, lcm_projection!, lcm_boundary
using Test
import Chevie

function d_fingerprints(P)
    out = Dict{Int, Tuple{Int, String, String}}()
    for n in 0:(Gauss._top_degree(P.tree) - 1)
        a = sort!(Gauss._differential_lines(P, n))
        isempty(a) && continue
        b = sort!(Gauss._differential_lines(P, n; bar = true))
        out[n] = (length(a), _sha_lines(a), _sha_lines(b))
    end
    return out
end

function reference_d_fingerprints(case, variant)
    rows, _ = read_fingerprints(joinpath(TREE_REFERENCE, "$(case)_$(variant).sha.csv"))
    return Dict(parse(Int, r.n) => (r.lines, r.atoms, r.bar) for r in rows if r.section == "D")
end

# ∂∂x = 0 for every cell x: ∂(λ σ [y] τ) = λ σ ∂(y) τ.
function square_zero(P)
    T = P.tree
    G = T.germ
    bd = [lcm_boundary(P, x) for x in 1:ncells(T)]
    for x in 1:ncells(T)
        acc = BimoduleTerm[]
        for t in bd[x], s in bd[t.cell]
            push!(acc, BimoduleTerm(s.cell, Gauss._mul(G, t.left, s.left), Gauss._mul(G, s.right, t.right), t.coef * s.coef))
        end
        isempty(Gauss._normalise!(acc)) || return false
    end
    return true
end

@testset "Projection onto the cells" begin
    @testset "B*(A₂) by hand" begin
        # Simples e, a₁, a₂, a₃, Δ = 0:4 with a₁a₂ = a₂a₃ = a₃a₁ = Δ; cells 1 = [],
        # 2 = [1], 3 = [2], 4 = [3], 5 = [2|3], 6 = [3|1].
        T = CellTree(dual_dihedral_germ(3))
        C = ChainTree(T)
        @test C.ncells == 6 && Gauss._nchains(C) == 6
        @test Gauss._child(C, 1, 1) == 2 && Gauss._child(C, 3, 3) == 5 && Gauss._child(C, 4, 1) == 6
        @test Gauss._child(C, 1, 4) == 0
        @test Gauss._child!(C, 1, 4) == 7 && Gauss._child!(C, 1, 4) == 7 && C.degree[7] == 1
        @test Gauss._letters(C, 5) == [2, 3] && Gauss._prefixes(C, 5) == [1, 3]
        @test_throws ArgumentError Gauss._child!(C, 1, 0)
        @test lcm_splitting(T, 1, 4) == (1, 2)          # [Δ]: redundant, matched with [1|2]
        @test lcm_splitting(T, 4, 1) == (1, 0)          # [3|1]: the cell
        @test lcm_splitting(T, 2, 2) == (0, 2)          # [1|2]: collapsible
        @test_throws ArgumentError lcm_splitting(T, 3, 1)   # a₂a₁ is not simple
        P = LCMProjectionData(T)
        # p([Δ]) = a₁·[2] + [1]·a₂, from the faces of [1|2] but the first.
        @test collect(lcm_projection!(P, [4])) == [BimoduleTerm(2, 0, 2, 1), BimoduleTerm(3, 1, 0, 1)]
        @test isempty(lcm_projection!(P, [1, 2]))
        @test collect(lcm_projection!(P, [3, 1])) == [BimoduleTerm(6, 0, 0, 1)]
        @test isempty(lcm_boundary(P, 1))
        @test_throws ArgumentError lcm_boundary(P, 7)
        # The D lines of the tiny reference snapshot dual_A2_a.
        @test sort!(vcat(Gauss._differential_lines(P, 0), Gauss._differential_lines(P, 1))) == [
            "D 0 1 +1 1 - e", "D 0 1 -1 e - 1", "D 0 2 +1 2 - e", "D 0 2 -1 e - 2",
            "D 0 3 +1 3 - e", "D 0 3 -1 e - 3",
            "D 1 2.1 +1 2 3 e", "D 1 2.1 +1 e 2 3", "D 1 2.1 -1 1 2 e", "D 1 2.1 -1 e 1 2",
            "D 1 3.1 +1 3 1 e", "D 1 3.1 +1 e 3 1", "D 1 3.1 -1 1 2 e", "D 1 3.1 -1 e 1 2"]
        @test square_zero(P)
        @test startswith(sprint(show, P), "LCMProjectionData: 6 cells, ")
        @test startswith(sprint(show, P.chains), "ChainTree: ") && endswith(sprint(show, P.chains), " chains, of which 6 cells")
    end

    @testset "B⁺(A₂) by hand" begin
        P = LCMProjectionData(CellTree(classical_dihedral_germ(3)))
        # The D lines of the tiny reference snapshot classic_A2_a.
        @test sort!(vcat(Gauss._differential_lines(P, 0), Gauss._differential_lines(P, 1))) == [
            "D 0 1 +1 1 - e", "D 0 1 -1 e - 1", "D 0 2 +1 2 - e", "D 0 2 -1 e - 2",
            "D 1 2.1 +1 2 1 2", "D 1 2.1 +1 2.1 2 e", "D 1 2.1 +1 e 2 1.2",
            "D 1 2.1 -1 1 2 1", "D 1 2.1 -1 1.2 1 e", "D 1 2.1 -1 e 1 2.1"]
        @test square_zero(P)
    end

    @testset "Explicit germs: ∂² = 0" begin
        for G in (classical_dihedral_germ(5), dual_dihedral_germ(6), signed_interval_germ(4, :D, (2, -4, -3, 1)))
            for (v, (S, _)) in TREE_VARIANTS
                @test square_zero(LCMProjectionData(CellTree(G, S)))
                @test is_square_zero(lcm_resolution(G, S))
            end
        end
    end

    A(n) = Chevie.coxgroup(:A, n)
    B(n) = Chevie.coxgroup(:B, n)
    cases = [
        "classic_A3" => Chevie.BraidMonoid(A(3)),
        "classic_B3" => Chevie.BraidMonoid(B(3)),
        "dual_A3" => Chevie.DualBraidMonoid(A(3); c = collect(1:3)),
        "dual_B3" => Chevie.DualBraidMonoid(B(3); c = collect(1:3)),
        "dual_D4" => Chevie.DualBraidMonoid(Chevie.coxgroup(:D, 4); c = collect(1:4)),
        "G4" => Chevie.DualBraidMonoid(Chevie.crg(4)),
        "G5" => Chevie.DualBraidMonoid(Chevie.crg(5)),
    ]
    @testset "D lines of J0: $case" for (case, M) in cases
        G = germ(M)
        classical = startswith(case, "classic_")
        for (v, (S, form)) in TREE_VARIANTS
            P = LCMProjectionData(CellTree(G, S))
            @test square_zero(P)
            # The witness of the check of the certificates, which reads the matrices.
            @test is_square_zero(lcm_resolution(P))
            # Gauss_a of the classical monoids and Gauss_b: ∂² = 0 only (decision of 23 October).
            (classical && v == "a" || form === :counts) && continue
            mine, ref = d_fingerprints(P), reference_d_fingerprints(case, v)
            @test sort(collect(keys(mine))) == sort(collect(keys(ref)))
            k = form === :atoms ? 2 : form === :bar ? 3 : 1
            @test all(haskey(mine, n) && mine[n][1] == ref[n][1] && mine[n][k] == ref[n][k] for n in keys(ref))
        end
    end

    @testset "Lazy germ" begin
        # B⁺(A₃) from the lazy germ: the same D lines as from the tabulated germ.
        L = LCMProjectionData(CellTree(Gauss.CoxeterGerm(Chevie.cartan(A(3)))))
        P = LCMProjectionData(CellTree(germ(Chevie.BraidMonoid(A(3)))))
        @test all(sort(Gauss._differential_lines(L, n)) == sort(Gauss._differential_lines(P, n)) for n in 0:2)
        @test square_zero(L)
    end

    @testset "Exports" begin
        for s in (:ChainTree, :BimoduleTerm, :LCMProjectionData, :lcm_splitting, :lcm_projection!, :lcm_boundary)
            @test Base.ispublic(Gauss, s) && !Base.isexported(Gauss, s)
        end
    end
end
