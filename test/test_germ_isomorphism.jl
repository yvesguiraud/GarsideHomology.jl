# Tests of germ_isomorphism (src/germs/isomorphism.jl): isomorphisms found and checked on
# germs whose atoms are permuted, non-isomorphic germs with the same numbers of simples
# and atoms, a lazy germ, the limit of the search. The germs of Chevie (G(m, 1, n) and
# Bₙ, the root datum of GL_n) are in test_chevie_germs.jl.

isdefined(Main, :check_germ) || include("explicit_germ.jl")

using GarsideHomology
using Test

# A pseudo-random permutation of 1:n, reproducible from its seed.
function _isomorphism_test_permutation(n, seed)
    x = UInt64(seed)
    values = [(x = x * 0x5851f42d4c957f2d + 0x14057b7ef767814f) for _ in 1:n]
    return sortperm(values)
end

# The germ T with its atoms permuted: the atom k of the result is the atom perm[k] of T.
permuted_atoms(T, perm) = GarsideHomology.TableGerm(right_covers_of(T)[perm, :], 1)

# Whether φ is an isomorphism of germs from G to H, checked independently of
# germ_isomorphism: the image of each simple of G is computed along its smallest word,
# and every product by an atom is compared.
function is_germ_isomorphism(G, H, φ)
    length(φ) == natoms(G) == natoms(H) || return false
    f = Vector{Int32}(undef, nsimples(G))
    for s in GarsideHomology.simples(G)
        x = GarsideHomology.one(H)
        for a in GarsideHomology.word(G, s)
            x = GarsideHomology.product(H, x, φ[a])
            x === nothing && return false
        end
        f[s + 1] = x
    end
    allunique(f) || return false
    for s in GarsideHomology.simples(G), a in 1:natoms(G)
        p = GarsideHomology.product(G, s, a)
        q = GarsideHomology.product(H, f[s + 1], φ[a])
        (p === nothing) == (q === nothing) || return false
        p === nothing || f[p + 1] == q || return false
    end
    return true
end

const ISOMORPHISM_TEST_GERMS = [
    "B⁺(A₂)" => classical_dihedral_germ(3), "B⁺(I₂(4))" => classical_dihedral_germ(4),
    "B*(I₂(5))" => dual_dihedral_germ(5), "B*(I₂(6))" => dual_dihedral_germ(6),
    "B*(B₃)" => signed_interval_germ(3, :B, (2, 3, -1)),
    "[1, -1] in B₃" => signed_interval_germ(3, :B, (-1, -2, -3)),
    "B*(D₄)" => signed_interval_germ(4, :D, (2, -4, -3, 1)),
]

@testset "germ_isomorphism" begin

    @testset "Permuted atoms: $name" for (name, E) in ISOMORPHISM_TEST_GERMS
        T = GarsideHomology.TableGerm(E)
        φ = GarsideHomology.germ_isomorphism(T, T)
        @test φ isa Vector{Int32} && is_germ_isomorphism(T, T, φ)
        @test GarsideHomology.germ_isomorphism(E, T) !== nothing      # an ExplicitGerm and its table
        for seed in 1:3
            perm = _isomorphism_test_permutation(natoms(T), seed)
            P = permuted_atoms(T, perm)
            φ = GarsideHomology.germ_isomorphism(P, T)
            @test φ !== nothing && is_germ_isomorphism(P, T, φ)
            ψ = GarsideHomology.germ_isomorphism(T, P)
            @test ψ !== nothing && is_germ_isomorphism(T, P, ψ)
        end
    end

    @testset "Non-isomorphic germs" begin
        B3 = GarsideHomology.TableGerm(signed_interval_germ(3, :B, (2, 3, -1)))
        I3 = GarsideHomology.TableGerm(signed_interval_germ(3, :B, (-1, -2, -3)))
        # Same numbers of simples (20) and of atoms (9): told apart by the search.
        @test nsimples(B3) == nsimples(I3) && natoms(B3) == natoms(I3)
        @test GarsideHomology.germ_isomorphism(B3, I3) === nothing
        @test GarsideHomology.germ_isomorphism(I3, B3) === nothing
        # Different numbers of simples or of atoms: told apart without search.
        @test GarsideHomology.germ_isomorphism(GarsideHomology.TableGerm(classical_dihedral_germ(4)),
                                     GarsideHomology.TableGerm(dual_dihedral_germ(4))) === nothing
        @test GarsideHomology.germ_isomorphism(GarsideHomology.TableGerm(dual_dihedral_germ(5)),
                                     GarsideHomology.TableGerm(dual_dihedral_germ(6))) === nothing
        # Same counts and isomorphic: B⁺(I₂(2)) and B*(I₂(2)), two germs of A₁ × A₁, with
        # 4 simples and 2 atoms.
        @test GarsideHomology.germ_isomorphism(GarsideHomology.TableGerm(classical_dihedral_germ(2)),
                                     GarsideHomology.TableGerm(dual_dihedral_germ(2))) !== nothing
    end

    @testset "Lazy germ" begin
        L = GarsideHomology.CoxeterGerm([2 -1; -1 2])                 # B⁺(A₂), numbered on demand
        T = GarsideHomology.TableGerm(classical_dihedral_germ(3))
        φ = GarsideHomology.germ_isomorphism(L, T)
        @test φ !== nothing && GarsideHomology.nnumbered(L) == nsimples(L) == 6
        @test is_germ_isomorphism(T, L, GarsideHomology.germ_isomorphism(T, L))
    end

    @testset "Limit of the search" begin
        T = GarsideHomology.TableGerm(signed_interval_germ(3, :B, (2, 3, -1)))
        P = permuted_atoms(T, _isomorphism_test_permutation(natoms(T), 1))
        @test_throws "search limit 1 reached" GarsideHomology.germ_isomorphism(P, T; limit = 1)
        @test_throws ArgumentError GarsideHomology.germ_isomorphism(P, T; limit = 0)
        @test GarsideHomology.germ_isomorphism(P, T; limit = 1000) !== nothing
    end
end
