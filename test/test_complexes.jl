# Tests of the complexes layer (Gauss.Complexes): monoid algebras, based free modules,
# labelled sparse matrices, chain complexes, finitely generated modules. Small complexes
# written by hand, no monoid: a monoid here is only a list of names, its labels are
# integers.

using Gauss
using Gauss: rank     # public, not exported (Chevie exports its own `rank`)
using Test

# AbstractAlgebra and SparseArrays are dependencies of Gauss, not of the test
# environment: reach them through the submodule.
const AA = Gauss.Complexes.AbstractAlgebra
const SA = Gauss.Complexes.SparseArrays
const ZZ = AA.ZZ
const L2 = NTuple{2, Int32}

# A monoid described by the names of its simples; label l names simples[l + 1], so that
# label 0 is the identity.
struct NamedMonoid <: Gauss.Complexes.AbstractMonoid
    simples::Vector{String}
end

# A chain complex of another type, known to the layer only through the interface of
# AbstractChainComplex.
struct WrappedComplex{R, D, C} <: Gauss.Complexes.AbstractChainComplex{R, D}
    inner::C
end
WrappedComplex(C::ChainComplex{R, M, D}) where {R, M, D} = WrappedComplex{R, D, typeof(C)}(C)
Base.range(W::WrappedComplex) = range(W.inner)
Base.getindex(W::WrappedComplex, n::Integer) = W.inner[n]
Gauss.Complexes.differential(W::WrappedComplex, n::Integer) = differential(W.inner, n)
AA.base_ring(W::WrappedComplex) = base_ring(W.inner)

@testset "Complexes" begin

    @testset "MonoidAlgebra" begin
        m = NamedMonoid(["e", "a"])
        A = MonoidAlgebra(ZZ, m)
        @test base_ring(A) === ZZ
        @test A.side === :both
        @test Gauss.Complexes.label_type(A) == L2
        @test Gauss.Complexes.label_type(MonoidAlgebra(ZZ, m, :right)) == Int32
        @test_throws ArgumentError MonoidAlgebra(ZZ, m, :middle)
        @test_throws MethodError MonoidAlgebra(ZZ, :not_a_monoid)
        @test_throws MethodError MonoidAlgebra(1, m)
        @test A == MonoidAlgebra(ZZ, m, :both)
        @test hash(A) == hash(MonoidAlgebra(ZZ, m, :both))
        @test A != MonoidAlgebra(ZZ, m, :left)
        @test A != MonoidAlgebra(ZZ, NamedMonoid(["e", "b"]))
        @test Gauss.Complexes.identity_label(Int32) === Int32(0)
        @test Gauss.Complexes.identity_label(L2) === (Int32(0), Int32(0))
        @test occursin("bimodules", sprint(show, A))
    end

    @testset "BasedFreeModule" begin
        M = BasedFreeModule(ZZ, ["[]", "[a]"])
        @test rank(M) == 2
        @test base_ring(M) === ZZ
        @test cells(M) == ["[]", "[a]"]
        @test M == BasedFreeModule(ZZ, ["[]", "[a]"])
        @test hash(M) == hash(BasedFreeModule(ZZ, ["[]", "[a]"]))
        @test M != BasedFreeModule(ZZ, ["[a]", "[]"])
        N = BasedFreeModule(ZZ, 3)
        @test rank(N) == 3
        @test cells(N) == 1:3
        @test rank(BasedFreeModule(ZZ, 0)) == 0
        @test_throws ArgumentError BasedFreeModule(ZZ, -1)
        @test_throws MethodError BasedFreeModule(1, 3)
        @test N isa Gauss.Complexes.AbstractModule
        @test sprint(show, N) == "BasedFreeModule of rank 3 over $(ZZ)"
        A = MonoidAlgebra(ZZ, NamedMonoid(["e"]))
        @test base_ring(BasedFreeModule(A, 2)) === A
    end

    @testset "LabelledSparseMatrix" begin
        # Monomials (row, column, label, coefficient) in any order, with a duplicate to
        # sum (row 2, label 3) and a pair that cancels (row 3, label 2).
        I = [2, 1, 2, 2, 3, 3, 1]
        J = [1, 1, 1, 1, 2, 2, 2]
        L = Int32[3, 1, 1, 3, 2, 2, 2]
        V = [1, 2, -1, 4, 1, -1, 7]
        d = labelled_sparse(I, J, L, V, 3, 2)
        @test d == LabelledSparseMatrix(3, 2, [1, 4, 5], Int32[1, 2, 2, 1], Int32[1, 1, 3, 2],
                                        [2, -1, 5, 7])
        @test hash(d) == hash(LabelledSparseMatrix(3, 2, [1, 4, 5], Int32[1, 2, 2, 1],
                                                   Int32[1, 1, 3, 2], [2, -1, 5, 7]))
        @test d isa LabelledSparseMatrix{Int32, Int}
        @test size(d) == (3, 2) && size(d, 1) == 3 && size(d, 2) == 2 && size(d, 3) == 1
        @test nrows(d) == 3 && ncols(d) == 2
        @test SA.nnz(d) == 4
        @test SA.nzrange(d, 1) == 1:3 && SA.nzrange(d, 2) == 4:4
        @test SA.rowvals(d) == Int32[1, 2, 2, 1]
        @test labels(d) == Int32[1, 1, 3, 2]
        @test SA.nonzeros(d) == [2, -1, 5, 7]
        @test !iszero(d)
        @test occursin("3×2 LabelledSparseMatrix{Int32, Int64} with 4 stored entries",
                       sprint(show, d))

        z = LabelledSparseMatrix{Int32, Int}(3, 2)
        @test iszero(z) && size(z) == (3, 2) && SA.nnz(z) == 0
        @test labelled_sparse(Int[], Int[], Int32[], Int[], 3, 2) == z
        @test size(LabelledSparseMatrix{L2, Int}(0, 0)) == (0, 0)
        @test_throws ArgumentError LabelledSparseMatrix{Int32, Int}(-1, 2)

        # Pairs of labels are sorted lexicographically.
        e = labelled_sparse([1, 1], [1, 1], L2[(2, 1), (1, 2)], [1, -1], 1, 1)
        @test labels(e) == L2[(1, 2), (2, 1)] && SA.nonzeros(e) == [-1, 1]

        # Invariants.
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [1, 2], Int32[2], Int32[1], [1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [1, 2], Int32[1], Int32[-1], [1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [1, 2], Int32[1], Int32[1], [0])
        @test_throws ArgumentError LabelledSparseMatrix(2, 1, [1, 3], Int32[2, 1], Int32[1, 1], [1, 1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [1, 3], Int32[1, 1], Int32[1, 1], [1, 1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 2, [1, 2], Int32[1], Int32[1], [1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [2, 2], Int32[1], Int32[1], [1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 2, [1, 3, 2], Int32[1], Int32[1], [1])
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [1, 2], Int32[1], Int32[1], [1, 2])
        @test_throws ArgumentError LabelledSparseMatrix(1, 1, [1, 2], Int32[1], L2[(1, -1)], [1])
        # Label 0, the identity, is valid.
        @test labels(LabelledSparseMatrix(1, 1, [1, 2], Int32[1], L2[(0, 0)], [1])) == L2[(0, 0)]
        bad = LabelledSparseMatrix(1, 1, [1, 2], Int32[1], Int32[1], [0]; check = false)
        @test bad isa LabelledSparseMatrix
        @test_throws ArgumentError Gauss.Complexes.check_invariants(bad)
        @test Gauss.Complexes.check_invariants(d) === d

        # Monomials out of bounds, invalid labels, overflow.
        @test_throws ArgumentError labelled_sparse([2], [1], Int32[1], [1], 1, 1)
        @test_throws ArgumentError labelled_sparse([1], [2], Int32[1], [1], 1, 1)
        @test_throws ArgumentError labelled_sparse([1], [1], Int32[-1], [1], 1, 1)
        @test_throws ArgumentError labelled_sparse([1, 1], [1], Int32[1], [1], 1, 1)
        @test_throws OverflowError labelled_sparse([1, 1], [1, 1], Int32[1, 1],
                                                   [typemax(Int), 1], 1, 1)
        bd = labelled_sparse([1, 1], [1, 1], Int32[1, 1], BigInt[typemax(Int), 1], 1, 1)
        @test SA.nonzeros(bd) == [BigInt(typemax(Int)) + 1]
        # Coefficients other than machine integers: the arithmetic of their type.
        @test specialise(bd, l -> -1) == SA.sparse([1], [1], [-(BigInt(typemax(Int)) + 1)], 1, 1)
        @test !is_square_zero(ChainComplex([SA.sparse([1], [1], [BigInt(2)], 1, 1), SA.sparse([1], [1], [BigInt(3)], 1, 1)]))
    end

    @testset "Specialisation of a matrix" begin
        d = labelled_sparse([2, 1, 2, 2, 3, 3, 1], [1, 1, 1, 1, 2, 2, 2],
                            Int32[3, 1, 1, 3, 2, 2, 2], [1, 2, -1, 4, 1, -1, 7], 3, 2)
        A = specialise(d)
        @test A isa SA.SparseMatrixCSC{Int, Int}
        @test A == SA.sparse([1, 2, 1], [1, 1, 2], [2, 4, 7], 3, 2)
        χ(l) = isodd(l) ? 1 : -1
        @test specialise(d, χ) == SA.sparse([1, 2, 1], [1, 1, 2], [2, 4, -7], 3, 2)
        # Sums over labels that cancel leave no stored zero.
        e = labelled_sparse([1, 1], [1, 1], Int32[1, 2], [1, -1], 1, 1)
        @test SA.nnz(specialise(e)) == 0 && size(specialise(e)) == (1, 1)
        @test SA.nnz(specialise(e, χ)) == 1
        f = labelled_sparse([1, 1], [1, 1], Int32[1, 2], [typemax(Int), 1], 1, 1)
        @test_throws OverflowError specialise(f)
        @test size(specialise(LabelledSparseMatrix{Int32, Int}(0, 4))) == (0, 4)
    end

    @testset "Complexes over ℤ" begin
        # Circle: one vertex, one loop.
        circle = ChainComplex([SA.spzeros(Int, Int, 1, 1)])
        @test range(circle) == 0:1
        @test firstindex(circle) == 0 && lastindex(circle) == 1
        @test rank(circle[0]) == 1 && rank(circle[1]) == 1
        @test rank(circle[end]) == 1
        @test base_ring(circle) === ZZ
        @test differential(circle, 0) == SA.spzeros(Int, 1, 1)
        @test size(differential(circle, -1)) == (0, 1)
        @test size(differential(circle, 1)) == (1, 0)
        @test size(differential(circle, 5)) == (0, 0)
        @test differential(circle, 1) isa SA.SparseMatrixCSC{Int, Int}
        @test_throws BoundsError circle[2]
        @test is_square_zero(circle)
        @test sprint(show, circle) == "ChainComplex over $(ZZ) in degrees 0:1, ranks 1, 1"

        # Real projective plane: one cell in each degree, d₁ = 2.
        rp2 = ChainComplex([SA.spzeros(Int, Int, 1, 1), SA.sparse([1], [1], [2], 1, 1)])
        @test range(rp2) == 0:2
        @test differential(rp2, 1) == SA.sparse([1], [1], [2], 1, 1)

        # Torus: ranks 1, 2, 1, zero maps; shifted degrees.
        torus = ChainComplex([SA.spzeros(Int, Int, 1, 2), SA.spzeros(Int, Int, 2, 1)]; lo = 1)
        @test range(torus) == 1:3
        @test [rank(torus[n]) for n in range(torus)] == [1, 2, 1]
        @test size(differential(torus, 1)) == (1, 2) && size(differential(torus, 0)) == (0, 1)

        # 2-simplex: vertices 1, 2, 3, edges [12], [13], [23], one face.
        d0 = SA.sparse([1, 2, 1, 3, 2, 3], [1, 1, 2, 2, 3, 3], [-1, 1, -1, 1, -1, 1], 3, 3)
        d1 = SA.sparse([1, 2, 3], [1, 1, 1], [1, -1, 1], 3, 1)
        @test is_square_zero(ChainComplex([d0, d1]))
        @test !is_square_zero(ChainComplex([d0, SA.sparse([1, 2, 3], [1, 1, 1], [1, 1, 1], 3, 1)]))

        # Inconsistent data.
        @test_throws ArgumentError ChainComplex([SA.spzeros(Int, Int, 1, 2), SA.spzeros(Int, Int, 3, 1)])
        @test_throws ArgumentError ChainComplex(SA.SparseMatrixCSC{Int, Int}[])
        @test_throws ArgumentError ChainComplex(BasedFreeModule{typeof(ZZ), Base.OneTo{Int}}[],
                                                SA.SparseMatrixCSC{Int, Int}[])
        @test_throws ArgumentError ChainComplex([BasedFreeModule(ZZ, 1), BasedFreeModule(ZZ, 1)],
                                                SA.SparseMatrixCSC{Int, Int}[])
        @test_throws ArgumentError ChainComplex([BasedFreeModule(ZZ, 1), BasedFreeModule(ZZ, 2)],
                                                [SA.spzeros(Int, Int, 1, 1)])
        A = MonoidAlgebra(ZZ, NamedMonoid(["e"]))
        @test_throws ArgumentError ChainComplex([BasedFreeModule(ZZ, 1), BasedFreeModule(A, 1)],
                                                [LabelledSparseMatrix{L2, Int}(1, 1)])
        # A single module, no map.
        point = ChainComplex([BasedFreeModule(ZZ, 1)], SA.SparseMatrixCSC{Int, Int}[])
        @test range(point) == 0:0 && size(differential(point, 0)) == (1, 0)
    end

    @testset "Complexes over a monoid algebra" begin
        # Free monoid on one atom a, simples e = 0 and a = 1: resolution of bimodules
        # 0 → C₁ = ⟨[a]⟩ → C₀ = ⟨[]⟩, ∂[a] = a·[]·e − e·[]·a.
        A = MonoidAlgebra(ZZ, NamedMonoid(["e", "a"]))
        len = [0, 1]                                 # len[l + 1] = ℓ(label l)
        da = labelled_sparse([1, 1], [1, 1], L2[(1, 0), (0, 1)], [1, -1], 1, 1)
        P = ChainComplex([BasedFreeModule(A, ["[]"]), BasedFreeModule(A, ["[a]"])], [da])
        @test base_ring(P) == A
        @test differential(P, 0) === da
        @test differential(P, 1) == LabelledSparseMatrix{L2, Int}(1, 0)
        @test_throws ArgumentError is_square_zero(P)
        S = specialise(P)
        @test base_ring(S) === ZZ
        @test cells(S[1]) === cells(P[1])            # bases are shared, not copied
        @test iszero(differential(S, 0)) && size(differential(S, 0)) == (1, 1)
        # Sign on the right, trivial on the left: H₀(ℕ ; ℤε) = ℤ/2ℤ.
        sign_right(l) = isodd(len[l[2] + 1]) ? -1 : 1
        @test differential(specialise(P, sign_right), 0) == SA.sparse([1], [1], [2], 1, 1)

        # Free commutative monoid on a, b, simples e = 0, a = 1, b = 2, ab = 3: Koszul
        # resolution of bimodules, ∂[a|b] = (a⊗1 − 1⊗a)[b] − (b⊗1 − 1⊗b)[a].
        B = MonoidAlgebra(ZZ, NamedMonoid(["e", "a", "b", "ab"]))
        lenB = [0, 1, 1, 2]
        d0 = labelled_sparse([1, 1, 1, 1], [1, 1, 2, 2],
                             L2[(1, 0), (0, 1), (2, 0), (0, 2)], [1, -1, 1, -1], 1, 2)
        d1 = labelled_sparse([1, 1, 2, 2], [1, 1, 1, 1],
                             L2[(2, 0), (0, 2), (1, 0), (0, 1)], [-1, 1, 1, -1], 2, 1)
        T = ChainComplex([BasedFreeModule(B, ["[]"]), BasedFreeModule(B, ["[a]", "[b]"]),
                          BasedFreeModule(B, ["[a|b]"])], [d0, d1])
        @test sprint(show, T) == "ChainComplex over $(B) in degrees 0:2, ranks 1, 2, 1"
        ST = specialise(T)
        @test all(iszero(differential(ST, n)) for n in -1:2)
        @test is_square_zero(ST)
        sign_right_B(l) = isodd(lenB[l[2] + 1]) ? -1 : 1
        SE = specialise(T, sign_right_B)
        @test differential(SE, 0) == SA.sparse([1, 1], [1, 2], [2, 2], 1, 2)
        @test differential(SE, 1) == SA.sparse([1, 2], [1, 1], [-2, 2], 2, 1)
        @test is_square_zero(SE)
    end

    @testset "AbstractChainComplex" begin
        # Over ℤ: real projective plane and 2-simplex, wrapped.
        rp2 = ChainComplex([SA.spzeros(Int, Int, 1, 1), SA.sparse([1], [1], [2], 1, 1)])
        W = WrappedComplex(rp2)
        @test firstindex(W) == 0 && lastindex(W) == 2 && rank(W[end]) == 1
        @test homology(W) == homology(rp2) && homology(W, 1) == homology(rp2, 1)
        @test is_square_zero(W) && collapsible_pairs(W) == collapsible_pairs(rp2)
        d0 = SA.sparse([1, 2, 1, 3, 2, 3], [1, 1, 2, 2, 3, 3], [-1, 1, -1, 1, -1, 1], 3, 3)
        @test !is_square_zero(WrappedComplex(ChainComplex([d0, SA.sparse([1, 2, 3], [1, 1, 1], [1, 1, 1], 3, 1)])))
        # Over a monoid algebra: the Koszul resolution of the free commutative monoid on a, b.
        B = MonoidAlgebra(ZZ, NamedMonoid(["e", "a", "b", "ab"]))
        e0 = labelled_sparse([1, 1, 1, 1], [1, 1, 2, 2],
                             L2[(1, 0), (0, 1), (2, 0), (0, 2)], [1, -1, 1, -1], 1, 2)
        e1 = labelled_sparse([1, 1, 2, 2], [1, 1, 1, 1],
                             L2[(2, 0), (0, 2), (1, 0), (0, 1)], [-1, 1, 1, -1], 2, 1)
        T = ChainComplex([BasedFreeModule(B, ["[]"]), BasedFreeModule(B, ["[a]", "[b]"]),
                          BasedFreeModule(B, ["[a|b]"])], [e0, e1])
        V = WrappedComplex(T)
        @test V isa Gauss.Complexes.AbstractChainComplex{typeof(B), LabelledSparseMatrix{L2, Int}}
        sign_right(l) = isodd([0, 1, 1, 2][l[2] + 1]) ? -1 : 1
        for χ in (nothing, sign_right)
            S, ST = specialise(V, χ), specialise(T, χ)
            @test S isa ChainComplex && range(S) == range(ST)
            @test all(differential(S, n) == differential(ST, n) for n in -1:2)
            @test all(cells(S[n]) === cells(T[n]) for n in 0:2)
        end
        @test collapsible_pairs(V) == collapsible_pairs(T)
        @test_throws ArgumentError is_square_zero(V)
        @test_throws ArgumentError homology(V)
    end

    @testset "FinitelyGeneratedModule" begin
        G = FinitelyGeneratedAbelianGroup
        @test G === FinitelyGeneratedModule{BigInt}
        @test G(1) == G(1, Int[])
        @test torsion(G(0, [2, 3])) == [6]
        @test torsion(G(0, [6, 2])) == [2, 6]
        @test torsion(G(0, [4, 6])) == [2, 12]
        @test torsion(G(2, [1, 3, 1, 3, 6])) == [3, 3, 6]
        @test torsion(G(0, [BigInt(2)^70])) == [BigInt(2)^70]
        @test rank(G(2, [3])) == 2
        @test iszero(G(0)) && iszero(G(0, [1, 1]))
        @test !iszero(G(0, [2])) && !iszero(G(1))
        @test_throws ArgumentError G(-1)
        @test_throws ArgumentError G(0, [0])
        @test_throws ArgumentError G(0, [-2])
        @test_throws ArgumentError G(0, [2.0])
        @test_throws ArgumentError FinitelyGeneratedModule{Int}(1)   # only BigInt for now
        @test string(G(0)) == "0"
        @test string(G(1)) == "ℤ"
        @test string(G(2)) == "ℤ²"
        @test string(G(12)) == "ℤ¹²"
        @test string(G(0, [2])) == "ℤ/2ℤ"
        @test string(G(1, [3, 3, 6])) == "(ℤ/3ℤ)² ⊕ ℤ/6ℤ ⊕ ℤ"
        @test hash(G(1, [2])) == hash(parse(G, "ℤ/2ℤ ⊕ ℤ"))
        # Values in the forms of the plan (section 4.2) and of tests_homologie.txt.
        cases = [("0", 0, Int[]), ("ℤ", 1, Int[]), ("ℤ²", 2, Int[]), ("ℤ/2", 0, [2]),
                 ("ℤ/3 ⊕ ℤ", 1, [3]), ("ℤ/2 ⊕ ℤ/4", 0, [2, 4]), ("ℤ/2 ⊕ ℤ/6", 0, [2, 6]),
                 ("(ℤ/3)² ⊕ ℤ/6", 0, [3, 3, 6]), ("(ℤ/3)² ⊕ ℤ", 1, [3, 3]),
                 ("(ℤ/2)²", 0, [2, 2]), ("(ℤ/6)²", 0, [6, 6]), ("ℤ/3 ⊕ ℤ/6", 0, [3, 6]),
                 ("ℤ/3ℤ ⊕ ℤ/3ℤ ⊕ ℤ/6ℤ", 0, [3, 3, 6]), ("ℤ/2ℤ ⊕ ℤ", 1, [2]),
                 ("(ℤ/3ℤ)² ⊕ ℤ²", 2, [3, 3]), ("Z/2 + Z^2", 2, [2]), ("ℤ/2 ⊕ ℤ/3", 0, [6]),
                 ("ℤ/12ℤ ⊕ ℤ¹⁰", 10, [12])]
        for (s, r, t) in cases
            H = parse(G, s)
            @test H == G(r, t)
            @test parse(G, string(H)) == H
        end
        @test_throws ArgumentError parse(G, "ℤ/2ℤ ⊕ ℚ")
        @test_throws ArgumentError parse(G, "")
        @test_throws ArgumentError parse(G, "ℤ/2²")
    end

    @testset "Exports" begin
        for s in names(Gauss.Complexes)
            (s === :Complexes || !Base.isexported(Gauss.Complexes, s)) && continue
            @test Base.isexported(Gauss, s)
            @test getfield(Gauss, s) === getfield(Gauss.Complexes, s)
        end
        # Same functions as AbstractAlgebra's: no clash with AbstractAlgebra, Nemo, Hecke.
        @test nrows === AA.number_of_rows && ncols === AA.number_of_columns
        @test rank === AA.rank && base_ring === AA.base_ring
        # rank is public but not exported, because Chevie exports its own rank.
        @test !Base.isexported(Gauss, :rank) && Base.ispublic(Gauss, :rank)
        @test Gauss.rank === Gauss.Complexes.rank
        for s in (:AbstractMonoid, :AbstractModule, :AbstractChainComplex)
            @test Base.ispublic(Gauss.Complexes, s) && !Base.isexported(Gauss.Complexes, s)
        end
    end
end
