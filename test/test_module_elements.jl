# Tests of the elements of the free modules and of the differential as an object
# (src/complexes/elements.jl, J1.h), on modules written by hand: over ℤ with names as
# descriptors, over monoid algebras of a monoid known only by the numbers of its elements
# (the default hooks: labels written as numbers, no product but by the identity). No germ is
# built (rule of the layer of complexes).

using GarsideHomology
using GarsideHomology: MonoidAlgebra, BasedFreeModuleElement, Differential, rank
using Test

const _EL = GarsideHomology.Complexes
const _el_sparse = _EL.SparseArrays.sparse
const _EL_ZZ = _EL.AbstractAlgebra.ZZ

struct ElementsTestMonoid <: _EL.AbstractMonoid end
struct ElementsOtherMonoid <: _EL.AbstractMonoid end

@testset "Elements of free modules" begin
    M = BasedFreeModule(_EL_ZZ, ["a", "b", "c"])

    @testset "Over a ring" begin
        v = 2 * M("a") - M("c")
        @test v isa BasedFreeModuleElement && parent(v) === M && length(v) == 2 && !iszero(v)
        @test v == M([2, 0, -1]) == M("2[a] - [c]") == M("−[c] + 2[a]") && hash(v) == hash(M([2, 0, -1]))
        @test Vector(v) == [2, 0, -1] && _EL.SparseArrays.SparseVector(v) == _EL.SparseArrays.sparsevec([1, 3], [2, -1], 3)
        @test v["a"] == 2 && v["b"] == 0 && v["c"] == -1
        @test collect(v) == [(2, nothing, "a", nothing), (-1, nothing, "c", nothing)]
        @test string(v) == "2[a] - [c]" && string(zero(M)) == "0" && string(-M("b")) == "-[b]"
        @test iszero(zero(M)) && iszero(v - v) && iszero(0 * v) && v * 3 == 3 * v && zero(v) == zero(M)
        @test M("0") == zero(M) && M("b") == M("[b]") == M(" + [b] ") == M("1×[b]")
        @test M("2×[a]") == M("2*[a]") == M("2⋅[a]") == M("2[a]") == M("2([a])") == 2 * M("a")
        @test M("12[a] - 12[a] + [b]") == M("b")
        show_with(N, x) = sprint(show, x; context = :notation => N)
        @test show_with(Notation(coefficient_style = :times), v) == "2×[a] - [c]"
        @test show_with(Notation(coefficient_style = :star), v) == "2*[a] - [c]"
        @test show_with(Notation(coefficient_style = :parentheses), v) == "2([a]) - [c]"
        for bad in ("[d]", "d", "2[a]3", "x[a]", "[a", "2[a] -", "(2[a]", ")[a]", "(1)2[a]", "2([a])3", "2 a")
            @test_throws ArgumentError M(bad)
        end
        @test_throws ArgumentError M([1, 2])                       # not the rank
        @test_throws ArgumentError M(["a"])                        # not an element of the basis
        @test_throws ArgumentError M("a") + BasedFreeModule(_EL_ZZ, ["x"])("x")
        @test_throws OverflowError 2 * (typemax(Int) * M("a"))
    end

    @testset "Over monoid algebras" begin
        m = ElementsTestMonoid()
        B = BasedFreeModule(MonoidAlgebra(_EL_ZZ, m, :both), ["x", "y"])
        v = B("2×1[x]3 - [y]")
        @test string(v) == "2×1[x]3 - [y]" && v[1, "x", 3] == 2 && v[0, "y", 0] == -1 && v[0, "x", 0] == 0
        @test collect(v) == [(2, Int32(1), "x", Int32(3)), (-1, Int32(0), "y", Int32(0))]
        @test B("2[x]") == B("(2)[x]") != 2 * B("x") && string(B("2[x]")) == "2[x]"   # σ = 2, coefficient 1
        @test string(2 * B("[x]")) == "2×[x]" && B("2×[x]") == 2 * B("x") && B("(12)[x]") == B("12[x]")
        @test string(B("[x]1 + 1[x]")) == "1[x] + [x]1"                                # σ[x] before [x]τ
        @test_throws ArgumentError B("a[x]")                                          # not a number
        L = BasedFreeModule(MonoidAlgebra(_EL_ZZ, m, :left), ["x"])
        w = L("2×1[x]")
        @test string(w) == "2×1[x]" && w[1, "x"] == 2 && collect(w) == [(2, Int32(1), "x", nothing)]
        @test_throws ArgumentError L("[x]1")
        R = BasedFreeModule(MonoidAlgebra(_EL_ZZ, m, :right), ["x"])
        u = R("2[x]3")                                 # no σ on the right: 2 is the coefficient
        @test string(u) == "2[x]3" && u["x", 3] == 2 && collect(u) == [(2, nothing, "x", Int32(3))]
        @test R("1[x]") == R("x")                       # the coefficient 1, since σ cannot be there
        @test_throws ArgumentError R("(1)[x]")
        @test_throws ArgumentError _EL._act(v, :left, ElementsTestMonoid(), Int32(1))   # no product but by e
        @test _EL._act(B("[x]"), :left, m, Int32(1)) == B("1[x]")
        @test _EL._act(B("[x]"), :right, m, Int32(1)) == B("[x]1")
        @test _EL._act(w, :left, m, Int32(0)) == w
        @test_throws ArgumentError _EL._act(L("[x]"), :right, m, Int32(1))   # a left module
        @test_throws ArgumentError _EL._act(M("a"), :left, m, Int32(1))      # a ring
        @test_throws ArgumentError _EL._act(B("[x]"), :left, ElementsOtherMonoid(), Int32(1))   # another monoid
    end

    @testset "Differential" begin
        M0, M1 = BasedFreeModule(_EL_ZZ, ["p"]), BasedFreeModule(_EL_ZZ, ["a", "b"])
        C = ChainComplex([M0, M1], [_el_sparse([1 -1])])
        d = differential(C)
        @test d isa Differential && startswith(sprint(show, d), "Differential of ChainComplex")
        @test d(M1("a")) == M0("p") && d[M1("b")] == -M0("p") && iszero(d(M1("a") + M1("b")))
        @test iszero(d(M0("p"))) && rank(parent(d(M0("p")))) == 0      # the lowest degree
        @test iszero(d(d(M0("p")))) && iszero(d(BasedFreeModule(_EL_ZZ, 0)("0")))   # and below it
        @test_throws ArgumentError M("") 
        @test d(BasedFreeModule(_EL_ZZ, ["a", "b"])("a")) == M0("p")   # an equal module
        @test_throws ArgumentError d("a")                              # a descriptor has no degree
        @test_throws ArgumentError d(BasedFreeModule(_EL_ZZ, ["q"])("q"))
        twice = ChainComplex([M1, M1], [_el_sparse([1 0; 0 1])])
        @test_throws ArgumentError differential(twice)(M1("a"))       # in two degrees
        @test C("2[a] - [b]") == 2 * M1("a") - M1("b") && C("[p]") == M0("p")
        @test_throws ArgumentError C("[p] + [a]")                      # two degrees
        @test_throws ArgumentError C("[q]")
        # Over the bimodules: d[x] = 1[p] - [p]2; the boundary of σ[x] needs a product.
        m = ElementsTestMonoid()
        A = MonoidAlgebra(_EL_ZZ, m, :both)
        P0, P1 = BasedFreeModule(A, ["p"]), BasedFreeModule(A, ["x"])
        D = _EL.LabelledSparseMatrix([1, 1], [1, 1], [(Int32(1), Int32(0)), (Int32(0), Int32(2))], [1, -1], 1, 1)
        dB = differential(ChainComplex([P0, P1], [D]))
        @test string(dB(P1("x"))) == "1[p] - [p]2" && dB(P1("[x]")) == P0("1[p] - [p]2")
        @test_throws ArgumentError dB(P1("1[x]"))                     # 1·1 is not available
        # One side: σ·σ′ on the left, τ′·τ on the right.
        for (side, text) in ((:left, "1[p]"), (:right, "[p]1"))
            S = MonoidAlgebra(_EL_ZZ, m, side)
            Q0, Q1 = BasedFreeModule(S, ["p"]), BasedFreeModule(S, ["x"])
            E = _EL.LabelledSparseMatrix([1], [1], [Int32(1)], [1], 1, 1)
            @test differential(ChainComplex([Q0, Q1], [E]))(Q1("x")) == Q0(text)
        end
    end

    @testset "Exports" begin
        for s in (:BasedFreeModuleElement, :Differential)
            @test Base.ispublic(GarsideHomology, s) && !Base.isexported(GarsideHomology, s)
            @test Base.ispublic(_EL, s)
        end
    end
end
