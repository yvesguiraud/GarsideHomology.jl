# Tests of the elements of the modules of the resolutions and of their differential (J1.h):
# cells as generators, the action of the elements of the germ, writing and reading with the
# germ, d(x) against the columns of the matrices, d(d(x)) = 0, and the tensor products.

using GarsideHomology
using GarsideHomology: basis, rank, tensor_product, BasedFreeModuleElement, Differential
using Test
import Chevie

const _RE_SA = GarsideHomology.Complexes.SparseArrays

# Whether d(x) is, for every cell x of the complex C, the column of x in the matrix of d.
function _boundaries_are_columns(C)
    d = differential(C)
    for n in first(range(C)):(last(range(C)) - 1), (j, x) in enumerate(basis(C[n + 1]))
        D = differential(C, n)
        v = d(x)
        p = _RE_SA.nzrange(D, j)
        mine = sort(collect(zip(v.indices, v.labels, v.coefficients)))
        labs = D isa _RE_SA.SparseMatrixCSC ? fill(nothing, length(p)) : GarsideHomology.labels(D)[p]
        theirs = sort(collect(zip(_RE_SA.rowvals(D)[p], labs, _RE_SA.nonzeros(D)[p])))
        mine == theirs || return false
    end
    return true
end

@testset "Elements of resolutions" begin
    A2 = germ(Chevie.coxgroup(:A, 2))                    # B⁺(A₂): cells e; 1, 2; 2|12
    P = lcm_resolution(A2)
    d = differential(P)
    x = P("2|12")
    C2 = GarsideHomology.CoxeterGerm([2 -1; -1 2])       # the same monoid, another germ

    @testset "Generators and their boundaries" begin
        @test P(x) == P[2](x) && string(P(x)) == "[2|12]" && string(P(P())) == "[]"
        @test d isa Differential && d(P("1")) == d[P("1")] == d(P(P("1")))
        @test string(d(P("1"))) == "1[] - []1" && string(d(P("2"))) == "2[] - []2"
        @test iszero(d(d(x))) && parent(d(d(x))) == P[0]
        @test iszero(d(P())) && rank(parent(d(P()))) == 0
        @test P(string(d(x))) == d(x) && P[1](string(d(x))) == d(x)
        @test P("[]") == P("[e]") == P[0](P()) && P[0]("e") == P(P())
        @test_throws ArgumentError P[0]("")
        @test P("2×1[2]1 - [1]") == 2 * (A2(1) * P(P("2")) * A2(1)) - P(P("1"))
        @test string(P("2×1[2]1 - [1]")) == "-[1] + 2×1[2]1"
        @test_throws ArgumentError P("[1] + [2|12]")                 # two degrees
        @test_throws ArgumentError P[1]("[2|12]")                    # not of degree 1
        @test_throws ArgumentError P[1](x)
        @test_throws ArgumentError P("[1|2]")                        # not a cell
        v = d(P("1"))
        @test v[A2(1), P(), A2()] == 1 && v[A2(), P(), A2(1)] == -1 && v[A2(1), P(), A2(1)] == 0
        @test collect(v) == [(1, A2(1), P(), A2()), (-1, A2(), P(), A2(1))]
        for R in (P, lcm_resolution(germ(Chevie.coxgroup(:A, 3))), lcm_resolution(dual_germ(Chevie.coxgroup(:A, 3))),
                  lcm_resolution(germ(Chevie.coxgroup(:B, 3))), lcm_resolution(GarsideHomology._dual_dihedral_germ(65)))
            dR = differential(R)
            @test _boundaries_are_columns(R)
            @test all(iszero(dR(dR(c))) for n in range(R) for c in basis(R[n]))
            @test all(R(string(dR(c))) == dR(c) for n in range(R) for c in basis(R[n]) if !iszero(dR(c)))
        end
    end

    @testset "Actions of the germ" begin
        y = P(P("1"))
        @test A2(1) * P(P("2")) == P("1[2]") && P(P("1")) * A2(2) == P("[1]2")
        @test A2(2) * (A2(1) * y) == A2(2, 1) * y == P("21[1]")
        @test_throws ArgumentError A2(1) * (A2(1) * y)               # 1·1 is not a simple
        @test_throws ArgumentError C2(1) * y
        @test d(A2(2) * y) == A2(2) * d(y) == P("21[] - 2[]1")
        @test_throws ArgumentError y[C2(1), P(), A2()]
    end

    @testset "Settings of the display" begin
        try
            set_notation!(atom_names = ["a", "b"])
            v = 2 * d(P("2"))
            @test string(v) == "2b[] - 2[]b" && P("2b[] - 2[]b") == v
            @test string(d(x)) == string(d(P(x))) && P(string(d(x))) == d(x)
        finally
            set_notation!()
        end
        @test string(2 * d(P("2"))) == "2×2[] - 2×[]2" && P("2×2[] - 2×[]2") == 2 * d(P("2"))
        @test sprint(show, 2 * P(P("2")); context = :notation => Notation(coefficient_style = :parentheses)) == "2([2])"
        @test sprint(show, d(P("1")); context = :notation => Notation(empty_cell_in_module = "[e]")) == "1[e] - [e]1"
    end

    @testset "Tensor products" begin
        E = extend_scalars(P)
        dE = differential(E)
        @test E[1](P("1")) isa BasedFreeModuleElement && string(E(P("1"))) == "[1]"
        @test iszero(dE(P("1"))) && dE(x) == E[1](Vector(dE(x))) && _boundaries_are_columns(E)
        @test Vector(dE(x)) == Vector(differential(E, 1)[:, 1]) && _RE_SA.SparseVector(dE(x)) == differential(E, 1)[:, 1]
        @test string(2 * E(P("1"))) == "2[1]" && E("2[1] - [2]") == 2 * E(P("1")) - E(P("2"))
        R = tensor_product(ℤ, P)                                     # right modules: labels τ
        @test string(R("2[1]2")) == "2[1]2" && R("2[1]2") == 2 * (R(P("1")) * A2(2))
        @test_throws ArgumentError A2(1) * R(P("1"))
        @test _boundaries_are_columns(R)
        L = tensor_product(P, ℤ)                                     # left modules: labels σ
        @test string(L("2×1[2]")) == "2×1[2]" && L("2×1[2]") == 2 * (A2(1) * L(P("2")))
        @test_throws ArgumentError L(P("1")) * A2(1)
        @test _boundaries_are_columns(L)
        @test_throws ArgumentError dE(P(P("1")))                      # an element of P, not of E
    end
end
