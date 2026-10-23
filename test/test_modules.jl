# Tests of the modules of coefficients and of the tensor products (modules/): TrivialModule
# and SignModule, B*(A₂) and B*(G₄) by hand, one side and both sides, associativity,
# agreement with specialise, homology with coefficients, errors, lazy germ, notation.

isdefined(Main, :check_germ) || include("explicit_germ.jl")

using Gauss
using Gauss: MonoidAlgebra, AbstractMonoidModule, tensor_product, ⊗, cell_tree
using Test
import Chevie

const _M_AA = Gauss.Complexes.AbstractAlgebra
const _M_ZZ = _M_AA.ZZ
const _M_SA = Gauss.Complexes.SparseArrays

@testset "Modules of coefficients" begin
    G = FinitelyGeneratedAbelianGroup
    GG = FinitelyGeneratedGradedAbelianGroup
    trivial3 = GG([G(1), G(1), G(0)])                 # H_*(B₃; ℤ)
    sign3 = GG([G(0, [2]), G(0, [3]), G(0)])          # H_*(B₃; ℤε)
    T, S = TrivialModule(_M_ZZ), SignModule(_M_ZZ)
    R0 = lcm_resolution(dual_dihedral_germ(3))

    @testset "TrivialModule and SignModule" begin
        @test T isa AbstractMonoidModule{typeof(_M_ZZ)} && S isa AbstractMonoidModule{typeof(_M_ZZ)}
        @test base_ring(T) === _M_ZZ && Gauss.rank(T) == 1 && Gauss.rank(S) == 1
        @test T == TrivialModule(_M_ZZ) && S == SignModule(_M_ZZ) && T != S
        @test sprint(show, T) == "TrivialModule over $(_M_ZZ)" && sprint(show, S) == "SignModule over $(_M_ZZ)"
        @test base_ring(SignModule(_M_AA.QQ)) === _M_AA.QQ      # any ring, tensors over ℤ only
        @test_throws MethodError TrivialModule(2)
        # The integers by default: those of AbstractAlgebra, Gauss.ZZ (public, not exported).
        @test TrivialModule() == T && SignModule() == S && Gauss.ZZ === _M_ZZ
        @test Base.ispublic(Gauss, :ZZ) && !Base.isexported(Gauss, :ZZ)
        @test homology(R0, SignModule()) == sign3
    end

    R = lcm_resolution(dual_dihedral_germ(3))        # B*(A₂): cells [], [1], [2], [3], 2.1, 3.1

    @testset "B*(A₂) by hand" begin
        # One side: ℤ ⊗ R, labels τ; ∂[a] = [] − []·a, the terms of the same row and label summed.
        L = T ⊗ R
        @test L isa ChainComplex && base_ring(L) == MonoidAlgebra(_M_ZZ, cell_tree(R).germ, :right)
        @test cells(L[1]) === cells(R[1]) && cell_tree(L) === cell_tree(R)
        d0 = differential(L, 0)
        @test labels(d0) == Int32[0, 1, 0, 2, 0, 3] && _M_SA.nonzeros(d0) == [1, -1, 1, -1, 1, -1]
        d1 = differential(L, 1)
        @test _M_SA.rowvals(d1) == Int32[1, 2, 2, 3, 1, 1, 2, 3] && labels(d1) == Int32[2, 0, 3, 0, 0, 2, 0, 1]
        @test _M_SA.nonzeros(d1) == [-1, -1, 1, 1, 1, -1, -1, 1]
        @test base_ring(R ⊗ S) == MonoidAlgebra(_M_ZZ, cell_tree(R).germ, :left)
        # Both sides, ℤε on the right: d₀ = (2, 2, 2), d₁ = [1] − 2[2] + [3], 2[1] − [2] − [3].
        E = tensor_product(T, R, S)
        @test differential(E, 0) == _M_SA.sparse([1, 1, 1], [1, 2, 3], [2, 2, 2], 1, 3)
        @test differential(E, 1) == _M_SA.sparse([1, 2, 3, 1, 2, 3], [1, 1, 1, 2, 2, 2], [1, -2, 1, 2, -1, -1], 3, 2)
        @test homology(E) == sign3
        # Trivial coefficients: d₀ = 0, d₁ = [3] − [1], [3] − [2].
        Z = tensor_product(T, R, T)
        @test iszero(differential(Z, 0))
        @test differential(Z, 1) == _M_SA.sparse([1, 3, 2, 3], [1, 1, 2, 2], [-1, 1, -1, 1], 3, 2)
        @test homology(Z) == trivial3 == homology(R, T)
        # ℤε on the left only: another complex, the same homology; on both sides: trivial.
        @test homology(tensor_product(S, R, T)) == sign3
        @test homology(tensor_product(S, R, S)) == trivial3
        @test homology(R, S) == sign3 && homology(R, S, 1) == G(0, [3]) && homology(R, S, 0) == G(0, [2])
    end

    @testset "Associativity, agreement with specialise" begin
        for (N, M) in ((T, S), (S, T), (S, S), (T, T)), C in (R, lcm_resolution(dual_dihedral_germ(4), LCMOrdering(TrivialPreorder())))
            g = cell_tree(C).germ
            χ(m, s) = m isa SignModule && isodd(length(g, s)) ? -1 : 1
            direct = specialise(C, l -> χ(N, l[1]) * χ(M, l[2]))
            for X in ((N ⊗ C) ⊗ M, N ⊗ (C ⊗ M), tensor_product(N, C, M))
                @test X isa ChainComplex && base_ring(X) === _M_ZZ
                @test all(differential(X, n) == differential(direct, n) for n in range(C))
                @test all(cells(X[n]) === cells(C[n]) for n in range(C))
            end
        end
        # From one side: homology(ℤ ⊗ R, M) = H_*(B; M).
        @test homology(T ⊗ R, S) == sign3
    end

    @testset "B*(G₄)" begin
        # The dual germ of G₄ has 5 simples and 3 atoms, as that of B*(A₂); B(G₄) ≅ B₃.
        C = lcm_resolution(dual_germ(Chevie.crg(4)))
        @test natoms(cell_tree(C).germ) == 3 && nsimples(cell_tree(C).germ) == 5
        @test homology(C, T) == trivial3 && homology(C, S) == sign3
        # The germ of the monoid of Chevie (another numbering): the same homology.
        D = lcm_resolution(germ(Chevie.DualBraidMonoid(Chevie.crg(4))))
        @test homology(D, T) == trivial3 && homology(D, S) == sign3
    end

    @testset "Lazy germ" begin
        W = Chevie.coxgroup(:A, 3)
        lazy = lcm_resolution(Gauss.CoxeterGerm(Chevie.cartan(W)))
        table = lcm_resolution(germ(Chevie.BraidMonoid(W)))
        @test homology(lazy, T) == GG([G(1), G(1), G(0, [2]), G(0)]) == homology(table, T)   # H_*(B₄; ℤ)
        @test homology(lazy, S) == homology(table, S)
    end

    @testset "Errors" begin
        @test_throws ArgumentError T ⊗ (T ⊗ R)                     # left side already tensored
        @test_throws ArgumentError (R ⊗ S) ⊗ S                     # right side already tensored
        @test_throws ArgumentError tensor_product(T, T ⊗ R, S)
        @test_throws ArgumentError homology(R ⊗ S, S)              # complex of left modules
        @test_throws ArgumentError R ⊗ SignModule(_M_AA.QQ)        # not over the integers
        @test_throws ArgumentError specialise(R) ⊗ S              # not over a monoid algebra
        @test_throws ArgumentError S ⊗ specialise(R)
        @test_throws ArgumentError tensor_product(S, specialise(R), S)
        @test_throws ArgumentError homology(R)                     # still needs a module, or specialise
    end

    @testset "Exports" begin
        for s in (:TrivialModule, :SignModule)
            @test Base.isexported(Gauss, s)
        end
        for s in (:AbstractMonoidModule, :tensor_product, :⊗)
            @test Base.ispublic(Gauss, s) && !Base.isexported(Gauss, s)
        end
    end
end
