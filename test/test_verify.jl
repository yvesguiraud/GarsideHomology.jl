# Tests of verify and count_multichains (verify.jl), without Chevie: consistent germs,
# Fuss–Catalan numbers of dual germs, altered germs.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using Gauss
using Test

# A copy of a tabulated germ with one table replaced (the constructor of the structure
# does not check anything, which is what these tests need).
altered(T::Gauss.TableGerm{K}; lcms = T.atom_lcms, descents = T.descents) where {K} =
    Gauss.TableGerm{K}(T.natoms, T.lengths, descents, T.right_covers, T.left_quotients, lcms, T.lex_ranks)

@testset "verify and count_multichains" begin

    @testset "Consistent germs: $name" for (name, E) in [
            "B⁺(A₂)" => classical_dihedral_germ(3), "B⁺(I₂(6))" => classical_dihedral_germ(6),
            "B*(I₂(5))" => dual_dihedral_germ(5), "B*(I₂(70))" => dual_dihedral_germ(70),
            "[1, -1] in B₃" => signed_interval_germ(3, :B, (-1, -2, -3)),
            "B*(D₄)" => signed_interval_germ(4, :D, (2, -4, -3, 1))]
        T = Gauss.TableGerm(E)
        @test Gauss.verify(E) === nothing && Gauss.verify(E; full = true) === nothing
        @test Gauss.verify(T) === nothing && Gauss.verify(T; full = true) === nothing
        @test Gauss.count_multichains(T, 1) == nsimples(T) == Gauss.count_multichains(E, 1)
        @test Gauss.count_multichains(T, 2) == Gauss.count_multichains(E, 2)
    end

    @testset "Fuss–Catalan numbers" begin
        # Cat⁽ᵏ⁾(I₂(m)) = (2 + km)(m + km)/(2m) = (2 + km)(1 + k)/2.
        for m in (3, 5, 8), k in 1:3
            @test Gauss.count_multichains(Gauss.TableGerm(dual_dihedral_germ(m)), k) == (2 + k * m) * (1 + k) ÷ 2
        end
        # D₄: degrees 2, 4, 4, 6, h = 6.
        D4 = Gauss.TableGerm(signed_interval_germ(4, :D, (2, -4, -3, 1)))
        @test [Gauss.count_multichains(D4, k) for k in 1:3] == [50, 336, 1210]
        @test_throws ArgumentError Gauss.count_multichains(D4, 0)
    end

    @testset "Altered germs" begin
        T = Gauss.TableGerm(classical_dihedral_germ(4))      # e, s, t, st, ts, sts, tst, Δ = 0, …, 7
        @test Gauss.verify(T; full = true) === nothing
        # s ∨ t = Δ with a wrong complement: not a common multiple.
        lcms = copy(T.atom_lcms)
        lcms[2, 2] = (lcms[2, 2][1], Int32(0))
        @test_throws "not a common multiple" Gauss.verify(altered(T; lcms))
        # e ∨ s replaced by st: a common multiple, but not the least, seen only by the full check.
        st = Gauss.product(T, 1, 2)::Int32
        lcms = copy(T.atom_lcms)
        lcms[1, 1] = (st, st)
        @test Gauss.verify(altered(T; lcms)) === nothing
        @test_throws "does not divide the common multiple" Gauss.verify(altered(T; lcms); full = true)
        # A descent too many for st.
        descents = copy(T.descents)
        descents[st + 1] = Gauss.StaticBitSet{1}([1, 2])
        @test_throws "left_descents" Gauss.verify(altered(T; descents))
    end

    @test_throws MethodError Gauss.verify(Gauss.CoxeterGerm([2 -1; -1 2]))   # no simples
end
