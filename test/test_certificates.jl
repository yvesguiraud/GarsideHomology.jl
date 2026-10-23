# Tests of the certificates (certificates.jl) and of the Möbius functions of the germs
# (germs/verify.jl): B*(A₂) and B⁺(A₂) by hand, B*(B₃) in Gauss_d against the verdicts of J0
# (tasks C4 and C5), altered resolutions, the two computations of μ, the keywords, the lazy
# germ, the threads, is_square_zero.

isdefined(Main, :check_germ) || include("explicit_germ.jl")

using Gauss
using Gauss: LCMResolution, ResolutionCertificate, is_correct, is_minimal, is_linear,
             cell_tree, CellTree, CellBasis
using Test
import Chevie

const _CT_ZZ = Gauss.Complexes.AbstractAlgebra.ZZ
const _CT_SA = Gauss.Complexes.SparseArrays

# The fields of a certificate but the times.
_report(c::ResolutionCertificate) = [getfield(c, f) for f in fieldnames(ResolutionCertificate) if f !== :times]

# The resolution R with the matrix of d_n changed by f, which takes and returns the
# monomials (rows, columns, labels, values).
function _altered(f, R::LCMResolution, n::Integer)
    C = R.complex
    d = differential(R, n)
    I, J, L, V = Int[], Int[], NTuple{2, Int32}[], Int[]
    for j in 1:size(d, 2), k in _CT_SA.nzrange(d, j)
        push!(I, _CT_SA.rowvals(d)[k]); push!(J, j); push!(L, labels(d)[k]); push!(V, _CT_SA.nonzeros(d)[k])
    end
    I, J, L, V = f(I, J, L, V)
    maps = copy(C.maps)
    maps[n + 1] = labelled_sparse(I, J, L, V, size(d)...)
    return LCMResolution(cell_tree(R), ChainComplex(C.modules, maps))
end

# Drops the monomials number ks.
_drop(ks) = (I, J, L, V) -> (keep = setdiff(eachindex(I), ks); (I[keep], J[keep], L[keep], V[keep]))

@testset "Certificates" begin
    @testset "B*(A₂) by hand" begin
        G = dual_dihedral_germ(3)              # simples e, a₁, a₂, a₃, Δ = 0:4
        @test Gauss._mobius_values(G) == [1, -1, -1, -1, 2]
        @test Gauss._mobius_crosscut(G) == Dict{Int32, Int}(0 => 1, 1 => -1, 2 => -1, 3 => -1, 4 => 2)
        R = lcm_resolution(G)
        c = certify(R)
        @test c isa ResolutionCertificate && is_correct(c) && is_minimal(c) && is_linear(c)
        @test c.cells == [1, 3, 2] && c.terms == [6, 8]
        @test c.homogeneity == 0 && c.square_zero == 0 && c.euler == 0
        @test c.mobius_method === :crosscut && c.mobius_sum == 6 && c.split_weights == 0
        @test c.identity_terms == [0, 0] && c.max_label_length == [1, 1] && isempty(c.failures)
        @test first.(c.times) == [:homogeneity, :square_zero, :mobius, :euler]
        @test sprint(show, c) == "ResolutionCertificate: correct, minimal, linear"
        text = sprint(show, MIME"text/plain"(), c)
        @test startswith(text, "ResolutionCertificate: 6 cells, 14 terms\n  correct: yes")
        @test occursin("(μ by crosscut)", text) && occursin("minimal: yes (no term 1·y·1)", text)
        @test occursin("one degree per weight   yes (Σ|μ| = 6, 6 cells)", text)
        @test endswith(text, "linear: yes (max ℓ(σ) + ℓ(τ) by degree: 1, 1)")
        @test _report(certify(R; mobius = :layers))[[1:4; 6:end]] == _report(c)[[1:4; 6:end]]
        @test certify(R; mobius = :layers).mobius_method === :layers
        @test is_square_zero(R)
    end

    @testset "B⁺(A₂) by hand" begin
        G = classical_dihedral_germ(3)         # weak order of S₃
        expected = Dict(Int32[] => 1, Int32[1] => -1, Int32[2] => -1, Int32[1, 2] => 0, Int32[2, 1] => 0,
                        Int32[1, 2, 1] => 1)
        μ = Gauss._mobius_values(G)
        @test all(μ[s + 1] == expected[Gauss.word(G, s)] for s in simples(G))
        @test sort!(collect(values(Gauss._mobius_crosscut(G)))) == [-1, -1, 1, 1]
        c = certify(lcm_resolution(G))
        @test is_correct(c) && is_minimal(c) && !is_linear(c)
        @test c.cells == [1, 2, 1] && c.terms == [4, 6] && c.mobius_sum == 4 && c.split_weights == 0
        @test c.max_label_length == [1, 2]     # as in C5: classic_A2, 0:1;1:2
        @test sprint(show, c) == "ResolutionCertificate: correct, minimal, not linear (max ℓ(σ) + ℓ(τ) = 2)"
        @test endswith(sprint(show, MIME"text/plain"(), c), "linear: no (max ℓ(σ) + ℓ(τ) by degree: 1, 2)")
    end

    M = Chevie.DualBraidMonoid(Chevie.coxgroup(:B, 3); c = collect(1:3))
    GB = germ(M)

    @testset "B*(B₃) against J0" begin
        # Gauss_a: minimal and linear (C4: 38 cells, Σ|μ| = 38).
        a = certify(lcm_resolution(GB))
        @test is_correct(a) && is_minimal(a) && is_linear(a)
        @test sum(a.cells) == 38 && a.mobius_sum == 38 && sum(a.terms) == 162
        # Gauss_d: C4 gives 40 cells for Σ|μ| = 38, one cell in an unexpected degree; C5 gives
        # 167 terms, two of length 0, and the maxima 1, 2, 1.
        d = certify(lcm_resolution(GB, LCMOrdering(TrivialPreorder())))
        @test is_correct(d) && !is_minimal(d) && !is_linear(d)
        @test sum(d.cells) == 40 && d.mobius_sum == 38 && d.split_weights == 1
        @test sum(d.terms) == 167 && sum(d.identity_terms) == 2 && d.max_label_length == [1, 2, 1]
        @test sprint(show, d) == "ResolutionCertificate: correct, not minimal (2 terms 1·y·1), not linear (max ℓ(σ) + ℓ(τ) = 2)"
        text = sprint(show, MIME"text/plain"(), d)
        @test occursin("minimal: no (2 terms 1·y·1, by degree ", text)
        @test occursin("no, 1 weight in several degrees (Σ|μ| = 38, 40 cells)", text)
        @test !occursin("note:", text)
        @test Gauss._mobius_values(GB)[end] == -10 == Gauss._mobius_crosscut(GB)[Int32(nsimples(GB) - 1)]
    end

    @testset "Altered resolutions" begin
        G = dual_dihedral_germ(3)
        R = lcm_resolution(G)
        # A term removed (the first of the column of the cell 2.1 in d₁): ∂² ≠ 0, still homogeneous.
        c = certify(_altered(_drop([1]), R, 1))
        @test c.square_zero == 1 && c.homogeneity == 0 && c.euler == 0 && !is_correct(c)
        @test first.(c.failures) == [:square_zero]
        @test startswith(last(only(c.failures)), "d_0 ∘ d_1, cell 2.1: remaining term ")
        @test sprint(show, c) == "ResolutionCertificate: not correct (∂² ≠ 0 on 1 column), minimal, linear"
        @test occursin("\n      first: d_0 ∘ d_1, cell 2.1", sprint(show, MIME"text/plain"(), c))
        @test !is_square_zero(_altered(_drop([1]), R, 1))
        # A coefficient changed: ∂² ≠ 0.
        flip = (I, J, L, V) -> (W = copy(V); W[end] = -W[end]; (I, J, L, W))
        @test certify(_altered(flip, R, 1)).square_zero == 1
        # A label changed, a₁·[] into a₂·[] in the boundary of [a₁]: one inhomogeneous term.
        relabel = (I, J, L, V) -> (K = copy(L); k = findfirst(==((Int32(1), Int32(0))), K); K[k] = (2, 0); (I, J, K, V))
        h = certify(_altered(relabel, R, 0))
        @test h.homogeneity == 1 && !is_correct(h)
        @test Gauss._failure(h, :homogeneity) == "d_0, cell 1: term +1 2 [-] e of another weight"
        # The weight of a cell changed in a copy of the tree: the Euler characteristic fails on
        # its old and new weights.
        T = CellTree(G)
        T.weight[6] = T.weight[2]
        A = MonoidAlgebra(_CT_ZZ, G, :both)
        bad = LCMResolution(T, ChainComplex([BasedFreeModule(A, CellBasis(T, n)) for n in 0:2], R.complex.maps))
        e = certify(bad)
        @test e.euler == 2 && e.split_weights == 1 && e.homogeneity > 0 && !is_correct(e)
        @test occursin("inhomogeneous terms, Euler fails on 2 weights), minimal, linear", sprint(show, e))
        @test Gauss._failure(e, :euler) == "weight 1: Σₙ (−1)ⁿ cₙ = 0, μ(1, ω) = -1"
    end

    @testset "μ by layers and by crosscut" begin
        germs = [dual_dihedral_germ(5), classical_dihedral_germ(4), germ(Chevie.BraidMonoid(Chevie.coxgroup(:A, 3))),
                 germ(Chevie.DualBraidMonoid(Chevie.coxgroup(:A, 3); c = collect(1:3))), GB,
                 signed_interval_germ(4, :D, (2, -4, -3, 1)), germ(Chevie.DualBraidMonoid(Chevie.crg(4)))]
        for G in germs
            v = Gauss._mobius_values(G)
            d = Gauss._mobius_crosscut(G)
            @test all(v[s + 1] == get(d, s, 0) for s in simples(G)) && all(!iszero, values(d))
        end
        D4 = signed_interval_germ(4, :D, (2, -4, -3, 1))       # B*(D₄): μ(1, Δ) = 20
        lengths = [length(D4, s) for s in simples(D4)]
        @test Gauss._mobius_values(D4)[argmax(lengths)] == 20
    end

    @testset "Keywords" begin
        R = lcm_resolution(dual_dihedral_germ(3))
        @test_throws ArgumentError certify(R; mobius = :other)
        @test_throws ArgumentError certify(R; nworkers = 0)
        @test_throws ArgumentError certify(specialise(R))
        n = certify(R; mobius = :none)
        @test n.mobius_method === :none && n.euler === nothing && n.mobius_sum === nothing && n.split_weights === nothing
        @test !is_correct(n) && is_minimal(n) && first.(n.times) == [:homogeneity, :square_zero]
        @test sprint(show, n) == "ResolutionCertificate: not certified (Euler not checked), minimal, linear"
        @test occursin("one degree per weight   not checked", sprint(show, MIME"text/plain"(), n))
        s = certify(R; square_zero = false)
        @test s.square_zero === nothing && !is_correct(s) && s.euler == 0
        @test sprint(show, s) == "ResolutionCertificate: not certified (∂² not checked), minimal, linear"
        @test sprint(show, certify(R; square_zero = false, mobius = :none)) ==
              "ResolutionCertificate: not certified (∂² not checked, Euler not checked), minimal, linear"
        # 25 atoms: layers by default, crosscut refused.
        B5 = lcm_resolution(dual_germ(Chevie.coxgroup(:B, 5)))
        @test natoms(cell_tree(B5).germ) == 25
        cB5 = certify(B5)
        @test cB5.mobius_method === :layers && is_correct(cB5)
        @test_throws ArgumentError certify(B5; mobius = :crosscut)
    end

    @testset "Lazy germ" begin
        W = Chevie.coxgroup(:A, 3)
        L = lcm_resolution(Gauss.CoxeterGerm(Chevie.cartan(W)))
        c = certify(L)
        @test c.mobius_method === :crosscut && is_correct(c) && is_minimal(c) && !is_linear(c)
        @test sum(c.cells) == 8 && c.mobius_sum == 8 && sum(c.terms) == 36 && c.max_label_length == [1, 2, 4]
        @test _report(c) == _report(certify(lcm_resolution(germ(Chevie.BraidMonoid(W)))))
        @test_throws ArgumentError certify(L; mobius = :layers)
        @test is_square_zero(L)
    end

    @testset "Threads" begin
        R = lcm_resolution(GB, LCMOrdering(TrivialPreorder()))
        @test _report(certify(R; nworkers = 1)) == _report(certify(R; nworkers = 4))
        # Several columns of d₁ altered: same count and same first failure on 1 and 4 tasks.
        cols = (I, J, L, V) -> _drop(unique(k -> J[k], eachindex(J)))(I, J, L, V)
        X = _altered(cols, R, 1)
        one_task, four_tasks = certify(X; nworkers = 1), certify(X; nworkers = 4)
        @test one_task.square_zero == four_tasks.square_zero > 1
        @test one_task.failures == four_tasks.failures
        @test !is_square_zero(X; nworkers = 1) && !is_square_zero(X; nworkers = 4)
        @test is_square_zero(R; nworkers = 4)
    end

    @testset "Exports" begin
        @test Base.isexported(Gauss, :certify)
        for s in (:ResolutionCertificate, :is_correct, :is_minimal, :is_linear)
            @test Base.ispublic(Gauss, s) && !Base.isexported(Gauss, s)
        end
    end
end
