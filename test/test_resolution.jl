# Tests of the lcm resolution (resolution.jl) and of the writing of its snapshots
# (snapshots.jl): the type LCMResolution, structure of the complex and of its bases,
# homology of small cases, snapshots written, read back by the tool of J1.b and compared
# with the references.

isdefined(Main, :read_snapshot) || include("snapshot_io.jl")
isdefined(Main, :check_germ) || include("explicit_germ.jl")

using Gauss
using Gauss: CellTree, CellBasis, LCMResolution, cell_tree, LCMProjectionData
# write_snapshot is called as Gauss.write_snapshot: the tool of J1.b (snapshot_io.jl)
# defines a write_snapshot of its own in Main, for its Snapshot type.
using Test
import Chevie

@testset "lcm resolution" begin
    G = FinitelyGeneratedAbelianGroup
    GG = FinitelyGeneratedGradedAbelianGroup
    braid3 = GG([G(1), G(1), G(0)])                 # braid group on three strands
    artinB3 = GG([G(1), G(2), G(2), G(1)])          # Artin group of type B₃

    @testset "B*(A₂) and B⁺(A₂)" begin
        C = lcm_resolution(dual_dihedral_germ(3))
        @test C isa LCMResolution && C isa Gauss.Complexes.AbstractChainComplex
        @test range(C) == 0:2 && [rank(C[n]) for n in 0:2] == [1, 3, 2]
        @test firstindex(C) == 0 && lastindex(C) == 2 && rank(C[end]) == 2
        @test_throws BoundsError C[3]
        @test base_ring(C) == MonoidAlgebra(Gauss.Complexes.AbstractAlgebra.ZZ, cell_tree(C).germ, :both)
        @test cells(C[1]) isa CellBasis && collect(cells(C[2])) == [5, 6]
        @test sprint(show, cells(C[1])) == "CellBasis: 3 cells of degree 1"
        @test sprint(show, MIME"text/plain"(), cells(C[1])) == "CellBasis: 3 cells of degree 1"
        @test sprint(show, C) == "LCMResolution: 6 cells, by degree [1, 3, 2], 14 terms, " *
                                 "germ ExplicitGerm with 3 atoms, LCMOrdering(ComplementDeglex())"
        T = cell_tree(C)
        @test T isa CellTree && cell_tree(specialise(C)) === T
        # The functions of the complexes layer, through the interface of AbstractChainComplex.
        @test specialise(C) isa ChainComplex && cells(specialise(C)[1]) === cells(C[1])
        @test isempty(collapsible_pairs(C))
        @test_throws ArgumentError homology(C)
        @test_throws ArgumentError is_square_zero(C)
        @test differential(C, 5) == Gauss.Complexes.LabelledSparseMatrix{NTuple{2, Int32}, Int}(0, 0)
        # The constructor checks that the complex is that of the tree.
        @test_throws ArgumentError LCMResolution(cell_tree(lcm_resolution(dual_dihedral_germ(3))), C.complex)
        plain = ChainComplex([BasedFreeModule(base_ring(C), collect(cells(C[n]))) for n in 0:2], C.complex.maps)
        @test_throws ArgumentError LCMResolution(T, plain)
        left = ChainComplex([BasedFreeModule(MonoidAlgebra(Gauss.Complexes.AbstractAlgebra.ZZ, T.germ, :left), cells(C[n])) for n in 0:2], C.complex.maps)
        @test_throws ArgumentError LCMResolution(T, left)
        @test size(differential(C, 0)) == (1, 3) && size(differential(C, 1)) == (3, 2)
        @test homology(specialise(C)) == braid3
        @test homology(specialise(lcm_resolution(classical_dihedral_germ(3)))) == braid3
        # Dehornoy–Lafont for another order of the atoms: same homology.
        R = lcm_resolution(dual_dihedral_germ(3), LCMOrdering(TrivialPreorder(); root_order = [3, 2, 1]))
        @test homology(specialise(R)) == braid3
        # With data kept by the caller, whose memo stays.
        P = LCMProjectionData(CellTree(dual_dihedral_germ(3)))
        D = lcm_resolution(P)
        @test !isempty(P.arena) && all(differential(D, n) == differential(C, n) for n in 0:1)
        @test_throws ArgumentError cell_tree(ChainComplex([Gauss.Complexes.SparseArrays.sparse([1 1])]))
    end

    @testset "Snapshots of B*(B₃)" begin
        M = Chevie.DualBraidMonoid(Chevie.coxgroup(:B, 3); c = collect(1:3))
        C = lcm_resolution(germ(M))
        @test homology(specialise(C)) == artinB3
        dir = mktempdir()
        path = joinpath(dir, "dual_B3_a.txt.gz")
        @test Gauss.write_snapshot(path, C; case = "dual_B3", extra = ["libelle" => "B*(B₃)"]) == path
        s = read_snapshot(path)
        @test header_value(s, "format") == "cliche-gauss 1" && header_value(s, "cas") == "dual_B3"
        @test header_value(s, "ordre") == "LCMOrdering(ComplementDeglex())"
        @test header_value(s, "libelle") == "B*(B₃)" && header_value(s, "cellules") == "38"
        @test startswith(header_value(s, "empreinte_germe"), "GermFingerprint(9 atoms")
        # The fingerprints written, those recomputed from the file read back, and those of J0.
        rows = snapshot_fingerprints(s)
        written, file = read_fingerprints(joinpath(dir, "dual_B3_a.sha.csv"))
        reference, _ = read_fingerprints(joinpath(@__DIR__, "data", "reference", "dual_B3_a.sha.csv"))
        @test written == rows == reference
        @test file.size == filesize(path) && file.sha == bytes2hex(open(Gauss.SHA.sha256, path))
        # The same complex as the tiny reference snapshot of J0, up to the header.
        tiny = read_snapshot(joinpath(@__DIR__, "data", "snapshots", "dual_B3_a.txt.gz"))
        @test snapshot_fingerprints(tiny) == rows
        @test homology(specialise(s.complex)) == homology(specialise(tiny.complex)) == artinB3
        # Gauss_d, strict as well.
        Cd = lcm_resolution(germ(M), LCMOrdering(TrivialPreorder()))
        pd = joinpath(dir, "dual_B3_d.txt.gz")
        Gauss.write_snapshot(pd, Cd; fingerprints = false)
        @test !isfile(joinpath(dir, "dual_B3_d.sha.csv"))
        reference_d, _ = read_fingerprints(joinpath(@__DIR__, "data", "reference", "dual_B3_d.sha.csv"))
        @test snapshot_fingerprints(read_snapshot(pd)) == reference_d
        @test header_value(read_snapshot(pd), "cas") === nothing
        # The fingerprint of the germ given, or omitted.
        f = Gauss.fingerprint(germ(M))
        Gauss.write_snapshot(pd, Cd; fingerprints = false, germ_fingerprint = f)
        @test header_value(read_snapshot(pd), "empreinte_germe") == header_value(s, "empreinte_germe")
        Gauss.write_snapshot(pd, Cd; fingerprints = false, germ_fingerprint = nothing)
        @test header_value(read_snapshot(pd), "empreinte_germe") == "non calculée"
        Gauss.write_snapshot(pd, Cd; fingerprints = false, germ_fingerprint = "given by the caller")
        @test header_value(read_snapshot(pd), "empreinte_germe") == "given by the caller"
        @test_throws ArgumentError Gauss.write_snapshot(pd, Cd; germ_fingerprint = :other)
        # Errors.
        @test_throws ArgumentError Gauss.write_snapshot(joinpath(dir, "x.txt"), C)
        @test_throws ArgumentError Gauss.write_snapshot(joinpath(dir, "y.txt.gz"), specialise(C))
        @test_throws ArgumentError Gauss.write_snapshot(joinpath(dir, "z.txt.gz"), C; extra = ["a" => "b\nc"])
    end

    @testset "Example of the README" begin
        @test homology(specialise(lcm_resolution(dual_germ(Chevie.crg(4))))) == braid3
    end

    @testset "Lazy germ" begin
        C = lcm_resolution(Gauss.CoxeterGerm(Chevie.cartan(Chevie.coxgroup(:A, 3))))
        @test [rank(C[n]) for n in 0:3] == [1, 3, 3, 1]
        path = joinpath(mktempdir(), "classic_A3_a.txt.gz")
        Gauss.write_snapshot(path, C)
        s = read_snapshot(path)
        @test header_value(s, "empreinte_germe") == "non calculée (germe paresseux)"
        reference, _ = read_fingerprints(joinpath(@__DIR__, "data", "reference", "classic_A3_a.sha.csv"))
        @test snapshot_fingerprints(s) == reference
    end

    @testset "Exports" begin
        @test Base.isexported(Gauss, :lcm_resolution)
        for s in (:LCMResolution, :CellBasis, :cell_tree, :write_snapshot)
            @test Base.ispublic(Gauss, s) && !Base.isexported(Gauss, s)
        end
    end
end
