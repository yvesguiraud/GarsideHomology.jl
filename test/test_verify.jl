# Tests of verify, count_multichains and fingerprint (verify.jl), without Chevie:
# consistent germs, Fuss–Catalan numbers of dual germs, altered germs, fingerprints.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using GarsideHomology
using Test

# A copy of a tabulated germ with one table replaced (the constructor of the structure
# does not check anything, which is what these tests need).
altered(T::GarsideHomology.TableGerm{K}; lcms = T.atom_lcms, descents = T.descents, lengths = T.lengths) where {K} =
    GarsideHomology.TableGerm{K}(T.natoms, lengths, descents, T.right_covers, T.left_quotients, lcms, T.lex_ranks)

# The same germ with its atoms renumbered: atom a of the result is atom π[a] of T.
function permute_atoms(T, π)
    n, N = nsimples(T), natoms(T)
    R = zeros(Int32, N, n)
    for σ in 0:(n - 1), a in 1:N
        p = GarsideHomology.product(T, σ, π[a])
        R[a, σ + 1] = p === nothing ? 0 : p + 1
    end
    return GarsideHomology.TableGerm(R, 1)
end

@testset "verify, count_multichains and fingerprint" begin

    @testset "Consistent germs: $name" for (name, E) in [
            "B⁺(A₂)" => classical_dihedral_germ(3), "B⁺(I₂(6))" => classical_dihedral_germ(6),
            "B*(I₂(5))" => dual_dihedral_germ(5), "B*(I₂(70))" => dual_dihedral_germ(70),
            "[1, -1] in B₃" => signed_interval_germ(3, :B, (-1, -2, -3)),
            "B*(D₄)" => signed_interval_germ(4, :D, (2, -4, -3, 1))]
        T = GarsideHomology.TableGerm(E)
        @test GarsideHomology.verify(E) === nothing && GarsideHomology.verify(E; full = true) === nothing
        @test GarsideHomology.verify(T) === nothing && GarsideHomology.verify(T; full = true) === nothing
        @test GarsideHomology.count_multichains(T, 1) == nsimples(T) == GarsideHomology.count_multichains(E, 1)
        @test GarsideHomology.count_multichains(T, 2) == GarsideHomology.count_multichains(E, 2)
    end

    @testset "Fuss–Catalan numbers" begin
        # Cat⁽ᵏ⁾(I₂(m)) = (2 + km)(m + km)/(2m) = (2 + km)(1 + k)/2.
        for m in (3, 5, 8), k in 1:3
            @test GarsideHomology.count_multichains(GarsideHomology.TableGerm(dual_dihedral_germ(m)), k) == (2 + k * m) * (1 + k) ÷ 2
        end
        # D₄: degrees 2, 4, 4, 6, h = 6.
        D4 = GarsideHomology.TableGerm(signed_interval_germ(4, :D, (2, -4, -3, 1)))
        @test [GarsideHomology.count_multichains(D4, k) for k in 1:3] == [50, 336, 1210]
        @test_throws ArgumentError GarsideHomology.count_multichains(D4, 0)
    end

    @testset "Altered germs" begin
        T = GarsideHomology.TableGerm(classical_dihedral_germ(4))      # e, s, t, st, ts, sts, tst, Δ = 0, …, 7
        @test GarsideHomology.verify(T; full = true) === nothing
        # s ∨ t = Δ with a wrong complement: not a common multiple.
        lcms = copy(T.atom_lcms)
        lcms[2, 2] = (lcms[2, 2][1], Int32(0))
        @test_throws "not a common multiple" GarsideHomology.verify(altered(T; lcms))
        # e ∨ s replaced by st: a common multiple, but not the least, seen only by the full check.
        st = GarsideHomology.product(T, 1, 2)::Int32
        lcms = copy(T.atom_lcms)
        lcms[1, 1] = (st, st)
        @test GarsideHomology.verify(altered(T; lcms)) === nothing
        @test_throws "does not divide the common multiple" GarsideHomology.verify(altered(T; lcms); full = true)
        # A descent too many for st.
        descents = copy(T.descents)
        descents[st + 1] = GarsideHomology.StaticBitSet{1}([1, 2])
        @test_throws "left_descents" GarsideHomology.verify(altered(T; descents))
    end

    @test_throws MethodError GarsideHomology.verify(GarsideHomology.CoxeterGerm([2 -1; -1 2]))   # no simples

    @testset "Fingerprints" begin
        # μ(1, Δ), recomputed by brute force: m − 1 for B*(I₂(m)), (−1)² for B⁺(I₂(m)).
        for m in (3, 5, 8)
            @test GarsideHomology.fingerprint(GarsideHomology.TableGerm(dual_dihedral_germ(m))).mobius == m - 1
            @test GarsideHomology.fingerprint(GarsideHomology.TableGerm(classical_dihedral_germ(m))).mobius == 1
        end
        D4 = signed_interval_germ(4, :D, (2, -4, -3, 1))
        f = GarsideHomology.fingerprint(GarsideHomology.TableGerm(D4))
        @test f.natoms == 12 && f.delta_length == 4 && f.simples_by_length == [1, 12, 24, 12, 1]
        @test f.multichains == (336, 1210) && f.mobius == 20 && sum(f.descent_profile) == 50
        @test GarsideHomology.fingerprint(D4) == f && hash(GarsideHomology.fingerprint(D4)) == hash(f)
        @test startswith(sprint(show, f), "GermFingerprint(12 atoms, ℓ(Δ) = 4, 50 simples")
        # A permutation of the atoms: other tables, same fingerprint.
        for (E, π) in [(dual_dihedral_germ(5), [2, 1, 3, 4, 5]), (D4, [3, 1, 2, 12, 5, 6, 7, 8, 9, 10, 11, 4])]
            T = GarsideHomology.TableGerm(E)
            P = permute_atoms(T, π)
            @test GarsideHomology.fingerprint(P) == GarsideHomology.fingerprint(T)
            E === D4 || @test P != T
        end
        # [1, −1] in B₃: 20 simples and 9 atoms, like B*(B₃), but μ(1, Δ) = −7 instead of −10.
        @test GarsideHomology.fingerprint(GarsideHomology.TableGerm(signed_interval_germ(3, :B, (-1, -2, -3)))).mobius == -7
        # Two simples of maximal length (tst given the length of Δ in B⁺(I₂(4))).
        T4 = GarsideHomology.TableGerm(classical_dihedral_germ(4))
        lengths = copy(T4.lengths)
        lengths[7] = lengths[8]
        @test_throws "simples of maximal length" GarsideHomology.fingerprint(altered(T4; lengths))
    end
end
