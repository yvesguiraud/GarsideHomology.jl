# Tests of the collapsible pairs (pairs.jl): complexes written by hand, and the three tiny
# reference snapshots against the counts recorded in J0 (data/reference/collapsible_pairs.csv).
# The other snapshots are checked in test_snapshots.jl when GAUSS_SNAPSHOT_DIR is set.

isdefined(Main, :read_snapshot) || include("snapshot_io.jl")

using Gauss
using Test

const _P_SA = Gauss.Complexes.SparseArrays
const PAIR_REFERENCE = joinpath(@__DIR__, "data", "reference", "collapsible_pairs.csv")

struct PairsTestMonoid <: Gauss.Complexes.AbstractMonoid end

@testset "Collapsible pairs" begin
    L2 = NTuple{2, Int32}
    A = MonoidAlgebra(Gauss.Complexes.AbstractAlgebra.ZZ, PairsTestMonoid())

    @testset "Complexes written by hand" begin
        # ⟨a⟩: ∂[a] = a·[]·e − e·[]·a has no monomial 1·y·1.
        da = labelled_sparse([1, 1], [1, 1], L2[(1, 0), (0, 1)], [1, -1], 1, 1)
        P = ChainComplex([BasedFreeModule(A, 1), BasedFreeModule(A, 1)], [da])
        @test isempty(collapsible_pairs(P))
        @test isempty(identity_monomials(P))

        # C₁ = ⟨x₁, x₂, x₃⟩ → C₀ = ⟨y₁, y₂⟩ with ∂x₁ = y₁, ∂x₂ = 2·y₁ and
        # ∂x₃ = −y₁ + (1 + σ)·y₂: only (x₁, y₁) and (x₃, y₁) are collapsible. Labels: e = 0,
        # σ = 1.
        d = labelled_sparse([1, 1, 1, 2, 2], [1, 2, 3, 3, 3],
                            L2[(0, 0), (0, 0), (0, 0), (0, 0), (1, 0)], [1, 2, -1, 1, 1], 2, 3)
        C = ChainComplex([BasedFreeModule(A, 2), BasedFreeModule(A, 3)], [d])
        pairs = collapsible_pairs(C)
        @test pairs == [(n = 0, source = 1, target = 1), (n = 0, source = 3, target = 1)]
        @test pair_counts(pairs) == (2, 2, 1)
        monomials = identity_monomials(C)
        @test pair_counts(monomials) == (3, 4, 2)
        @test count(m -> abs(m.coefficient) != 1, monomials) == 1

        # Over ℤ, after specialisation: entries ±1.
        @test collapsible_pairs(specialise(C)) ==
              [(n = 0, source = 1, target = 1), (n = 0, source = 3, target = 1)]
        @test isempty(collapsible_pairs(ChainComplex([_P_SA.spzeros(Int, Int, 1, 1)])))
        rp2 = ChainComplex([_P_SA.spzeros(Int, Int, 1, 1), _P_SA.sparse([1], [1], [2], 1, 1)])
        @test isempty(collapsible_pairs(rp2))
        d0 = _P_SA.sparse([1, 2, 1, 3, 2, 3], [1, 1, 2, 2, 3, 3], [-1, 1, -1, 1, -1, 1], 3, 3)
        d1 = _P_SA.sparse([1, 2, 3], [1, 1, 1], [1, -1, 1], 3, 1)
        simplex_pairs = collapsible_pairs(ChainComplex([d0, d1]))
        @test length(simplex_pairs) == 9
        @test first(simplex_pairs) == (n = 0, source = 1, target = 1)
        @test last(simplex_pairs) == (n = 1, source = 1, target = 3)
    end

    @testset "Tiny reference snapshots" begin
        reference = read_pair_counts(PAIR_REFERENCE)
        @test length(reference) == 344
        @test reference["dual_B3_a"].targets !== nothing      # recorded with its targets
        for name in ("classic_A2_a", "dual_A2_a", "dual_B3_a")
            s = read_snapshot(joinpath(@__DIR__, "data", "snapshots", name * ".txt.gz"))
            monomials = identity_monomials(s.complex)
            @test pair_counts_agree(pair_counts(monomials), reference[name])
            @test length(collapsible_pairs(s.complex)) ==
                  count(m -> abs(m.coefficient) == 1, monomials)
        end
    end
end
