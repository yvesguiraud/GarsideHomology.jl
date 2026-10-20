# Tests of the snapshot tool (snapshot_io.jl): a round trip on a snapshot written by hand,
# the three tiny reference snapshots of data/snapshots/ against their fingerprints, and,
# when the environment variable GAUSS_SNAPSHOT_DIR names a folder of snapshots (CatLab
# j0/ref, j0/out/gros), every snapshot of that folder against its reference fingerprints.

include("snapshot_io.jl")

using Gauss
using Gauss: rank
using Test

const SNAP_REFERENCE = joinpath(@__DIR__, "data", "reference")
const SNAP_TINY = joinpath(@__DIR__, "data", "snapshots")
const SNAP_DIR = get(ENV, "GAUSS_SNAPSHOT_DIR", "")

# Free commutative monoid on two atoms a = 1, b = 2: Koszul resolution of bimodules,
# ∂[a] = a·[]·e − e·[]·a, ∂[b] likewise, ∂[a|b] = (a⊗1 − 1⊗a)[b] − (b⊗1 − 1⊗b)[a].
const KOSZUL_SNAPSHOT = """
# format: cliche-gauss 1
# cas: test_koszul
# libelle: ℕ² (free commutative monoid on two atoms)
# cellules: 4
# termes: 8
C 0 - [] e
C 1 1 [1] 1
C 1 2 [2] 2
C 2 1.2 [1|2] 1.2
D 0 1 +1 1 - e
D 0 1 -1 e - 1
D 0 2 +1 2 - e
D 0 2 -1 e - 2
D 1 1.2 +1 1 2 e
D 1 1.2 +1 e 1 2
D 1 1.2 -1 2 1 e
D 1 1.2 -1 e 2 1
"""

function _write_gzip(path, text)
    stream = GzipCompressorStream(open(path, "w"))
    try
        write(stream, text)
    finally
        close(stream)
    end
    return path
end

# Read a snapshot given as text.
_read_text(text) = read_snapshot(_write_gzip(tempname() * ".txt.gz", text))

# Sign on both sides, (σ, τ) ↦ (-1)^(ℓ(σ) + ℓ(τ)), for the simples of a snapshot.
_sign(m::SnapshotMonoid) = l -> iseven(simple_length(m, l[1]) + simple_length(m, l[2])) ? 1 : -1

@testset "Snapshots" begin
    workdir = mktempdir()

    @testset "Round trip on a snapshot written by hand" begin
        source = _write_gzip(joinpath(workdir, "test_koszul_a.txt.gz"), KOSZUL_SNAPSHOT)
        s = read_snapshot(source)
        @test header_value(s, "cas") == "test_koszul"
        @test header_value(s, "absent") === nothing
        @test s.monoid.simples == ["e", "1", "2", "1.2"]
        @test [simple_length(s.monoid, i) for i in 1:4] == [0, 1, 1, 2]
        C = s.complex
        @test range(C) == 0:2
        @test [rank(C[n]) for n in range(C)] == [1, 2, 1]
        @test cells(C[2]) == [SnapshotCell("1.2", "[1|2]", 4)]
        L2 = NTuple{2, Int32}
        @test differential(C, 0) ==
              labelled_sparse([1, 1, 1, 1], [1, 1, 2, 2], L2[(2, 1), (1, 2), (3, 1), (1, 3)],
                              [1, -1, 1, -1], 1, 2)
        @test differential(C, 1) ==
              labelled_sparse([1, 1, 2, 2], [1, 1, 1, 1], L2[(3, 1), (1, 3), (2, 1), (1, 2)],
                              [-1, 1, 1, -1], 2, 1)

        # Writing regenerates the same text, hence the same snapshot.
        target = write_snapshot(joinpath(workdir, "copy.txt.gz"), s)
        @test gunzip_text(target) == KOSZUL_SNAPSHOT
        @test snapshots_equal(read_snapshot(target), s)

        # Fingerprints: the cell of degree 0 is the same in every snapshot.
        fp = snapshot_fingerprints(s)
        @test [(r.section, r.n, r.lines) for r in fp] ==
              [("C", "0", 1), ("C", "1", 2), ("C", "2", 1), ("D", "0", 4), ("D", "1", 4), ("total", "", 12)]
        @test fp[1].atoms == "29f9bb68ac2c7f6ea7e48352923a20a3767a6e76c7368cdd8b60382dbc017e97"
        @test fp[1].bar == "daf3904c702a460cb27c19f1c4c5abe0b6e7fa138ebff94515cec5570412127e"

        # Specialisations: trivial on both sides, then the sign on the right.
        @test all(iszero(differential(specialise(C), n)) for n in 0:1)
        sign_right(l) = isodd(simple_length(s.monoid, l[2])) ? -1 : 1
        SE = specialise(C, sign_right)
        @test differential(SE, 0) == Gauss.Complexes.SparseArrays.sparse([1, 1], [1, 2], [2, 2], 1, 2)
        @test is_square_zero(SE)
    end

    @testset "Malformed snapshots" begin
        bad(changes...) = _read_text(replace(KOSZUL_SNAPSHOT, changes...))
        @test_throws ArgumentError bad("cliche-gauss 1" => "cliche-gauss 2")       # format
        @test_throws ArgumentError bad("# termes: 8" => "# termes: 9")             # counts
        @test_throws ArgumentError bad("# cellules: 4" => "# cellules: 5")
        @test_throws ArgumentError bad("D 0 1 +1 1 - e" => "D 0 1 +0 1 - e")       # zero coefficient
        @test_throws ArgumentError bad("D 0 1 +1 1 - e" => "D 0 1 +1 1 3 e")       # unknown target
        @test_throws ArgumentError bad("D 0 1 +1 1 - e" => "D 0 3 +1 1 - e")       # unknown source
        @test_throws ArgumentError bad("D 0 1 +1 1 - e" => "D 0 1.2 +1 1 - e")     # source of degree 2
        @test_throws ArgumentError bad("D 0 1 +1 1 - e" => "D 0 1 +1 1 - e x")     # 8 fields
        @test_throws ArgumentError bad("C 1 2 [2] 2" => "C 1 1 [2] 2")             # repeated cell
        @test_throws ArgumentError bad("C 2 1.2 [1|2] 1.2\n" => "")               # d_1 without C_2
        @test_throws ArgumentError bad("D 0 1 +1 1 - e" => "X 0 1 +1 1 - e")       # unknown line type
        @test_throws ArgumentError bad("# termes: 8" => "# termes: 9",             # repeated term
                                       "D 0 1 -1 e - 1\n" => "D 0 1 -1 e - 1\nD 0 1 -1 e - 1\n")
        @test_throws ArgumentError bad("D 1 1.2 -1 e 2 1\n" =>                     # cell after terms
                                       "D 1 1.2 -1 e 2 1\nC 3 1.2.1 [1|2|1] 1.2\n")
        @test_throws ArgumentError bad("D 1 1.2 -1 e 2 1\n" =>                     # d_0 after d_1
                                       "D 1 1.2 -1 e 2 1\nD 0 2 -1 e - 2\n")
    end

    @testset "Tiny reference snapshots: $name" for name in ("classic_A2_a", "dual_A2_a", "dual_B3_a")
        path = joinpath(SNAP_TINY, name * ".txt.gz")
        rows, file = read_fingerprints(joinpath(SNAP_REFERENCE, name * ".sha.csv"))
        # The copy in the repository is the file the fingerprints were computed from.
        @test filesize(path) == file.size
        @test bytes2hex(open(SHA.sha256, path)) == file.sha
        s = read_snapshot(path)
        @test snapshot_fingerprints(s) == rows
        target = write_snapshot(joinpath(workdir, name * ".txt.gz"), s)
        @test gunzip_text(target) == gunzip_text(path)
        @test snapshots_equal(read_snapshot(target), s)
        @test is_square_zero(specialise(s.complex))
        @test is_square_zero(specialise(s.complex, _sign(s.monoid)))
    end

    if isempty(SNAP_DIR)
        @info "Reference snapshots of J0 not checked: set GAUSS_SNAPSHOT_DIR to a folder " *
              "of snapshots (.txt.gz) to check them against data/reference."
    else
        @testset "Snapshots of GAUSS_SNAPSHOT_DIR" begin
            files = sort!(filter(f -> endswith(f, ".txt.gz"), readdir(SNAP_DIR)))
            @test !isempty(files)
            results = NamedTuple{(:name, :lines, :seconds, :same_gzip), Tuple{String, Int, Float64, Bool}}[]
            println("Checking $(length(files)) snapshots of $SNAP_DIR")
            elapsed = @elapsed for f in files
                name = f[1:(end - 7)]
                reference = joinpath(SNAP_REFERENCE, name * ".sha.csv")
                if !isfile(reference)
                    @warn "No reference fingerprints for $f"
                    continue
                end
                @testset "$name" begin
                    t = @elapsed begin
                        rows, file = read_fingerprints(reference)
                        s = read_snapshot(joinpath(SNAP_DIR, f))
                        @test snapshot_fingerprints(s) == rows
                        target = write_snapshot(joinpath(workdir, f), s)
                        same_gzip = bytes2hex(open(SHA.sha256, target)) == file.sha
                        rm(target)
                    end
                    push!(results, (name = name, lines = rows[end].lines, seconds = t,
                                    same_gzip = same_gzip))
                    println("  ", rpad(name, 22), lpad(rows[end].lines, 9), " lines ",
                            lpad(round(t; digits = 2), 8), " s", same_gzip ? "" : "  (gzip bytes differ)")
                end
            end
            @info "$(length(results)) snapshots checked in $(round(elapsed; digits = 1)) s; " *
                  "gzip bytes identical for $(count(r -> r.same_gzip, results))."
        end
    end
end
