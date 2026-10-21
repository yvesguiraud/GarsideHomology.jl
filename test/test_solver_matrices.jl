# Comparison of the specialised differentials of B*(G34), B*(G35), B*(G36) in Gauss_a with
# the matrices of SparseElementaryDivisors (its test/data/G<p>_d<n>.mtx, exported in J0
# from matrices/G<p>.jld2 of the experiment repository), up to permutation of the bases.
# The bases of these matrices are given by their bar forms, written by the script
# j1/b5_bases.jl of the experiment repository. Runs when three environment variables are set:
#   GAUSS_SOLVER_MATRICES  the folder test/data of SparseElementaryDivisors
#   GAUSS_SOLVER_BASES     the folder j1/out/bases of the experiment repository
#   GAUSS_SNAPSHOT_DIR     a folder with the snapshots G34_a, G35_a, G36_a (.txt.gz)

isdefined(Main, :read_snapshot) || include("snapshot_io.jl")

using Gauss
using Test

const SOLVER_MATRICES = get(ENV, "GAUSS_SOLVER_MATRICES", "")
const SOLVER_BASES = get(ENV, "GAUSS_SOLVER_BASES", "")
const SOLVER_SNAPSHOTS = get(ENV, "GAUSS_SNAPSHOT_DIR", "")
const _read_matrix_market = Gauss.Complexes.SparseElementaryDivisors.read_matrix_market

"The bar forms of a file of bases, by degree: lines `C <degree> <bar form>`."
function read_bases(path::AbstractString)
    bases = Dict{Int, Vector{String}}()
    for l in eachline(path)
        (isempty(l) || startswith(l, "#")) && continue
        f = split(l, ' ')
        length(f) == 3 && f[1] == "C" || throw(ArgumentError("$path: unexpected line $l"))
        push!(get!(bases, parse(Int, f[2]), String[]), String(f[3]))
    end
    return bases
end

if isempty(SOLVER_MATRICES) || isempty(SOLVER_BASES) || isempty(SOLVER_SNAPSHOTS)
    @info "Matrices of the solver not compared: set GAUSS_SOLVER_MATRICES, GAUSS_SOLVER_BASES " *
          "and GAUSS_SNAPSHOT_DIR (see test_solver_matrices.jl)."
else
    @testset "Matrices of the solver, B*(G$p)" for p in (34, 35, 36)
        snapshot = joinpath(SOLVER_SNAPSHOTS, "G$(p)_a.txt.gz")
        bases_file = joinpath(SOLVER_BASES, "G$(p)_bases.txt")
        @test isfile(snapshot) && isfile(bases_file)
        if isfile(snapshot) && isfile(bases_file)
            S = specialise(read_snapshot(snapshot).complex)
            solver = read_bases(bases_file)
            @test sort!(collect(keys(solver))) == collect(range(S))
            # Position of each cell in our basis of degree n, from its bar form.
            position = Dict(n => Dict(x.bar => i for (i, x) in enumerate(cells(S[n]))) for n in range(S))
            for n in range(S)
                @test sort(solver[n]) == sort([x.bar for x in cells(S[n])])
            end
            for n in first(range(S)):(last(range(S)) - 1)
                A = _read_matrix_market(joinpath(SOLVER_MATRICES, "G$(p)_d$(n).mtx"))
                rows = [position[n][b] for b in solver[n]]
                columns = [position[n + 1][b] for b in solver[n + 1]]
                @test differential(S, n)[rows, columns] == A
            end
            println("  B*(G$p): bases and differentials d_0 to d_$(last(range(S)) - 1) compared")
        end
    end
end
