# benchmark/construction.jl — one round of the timing campaign of the construction of the
# lcm resolution (J1.d, task d7), in a fresh Julia session.
#
# Usage (normally run by run_construction.jl):
#     julia --project=benchmark --threads=auto benchmark/construction.jl <machine> <case> <path> <round>
# <path>: P1, through the monoid of Chevie, as in J0 (DualBraidMonoid or BraidMonoid,
# then germ(M)); P2, the package alone (dual_germ(W), or the lazy germ CoxeterGerm(cartan(W))
# of a classical braid monoid); in both, the group W comes from Chevie.
# Warm-up on a small germ of the same type (compilation out of the measures), then each
# step timed by @timed: group, monoid (P1 only), germ, tree (CellTree), differential
# (LCMProjectionData and lcm_resolution, the memo dropped at the end); after each step,
# Sys.maxrss and the live memory after a full GC. Then, out of the measures, the
# fingerprint of a tabulated germ, to compare the germs of P1 and P2.
# Output: one line in benchmark/out/construction_<machine>.csv.

using GarsideHomology
using GarsideHomology: CellTree, LCMProjectionData, ncells
import Chevie

const OUT_DIR = joinpath(@__DIR__, "out")

# The seven reference cases of J0, Gauss_a, with the constructions of CatLab/j1/common.jl;
# K: number of 64-bit words of the descent sets (atoms ≤ 64: 1).
const CASES = Dict(
    "G34" => (kind = :complex, K = 2, group = () -> Chevie.crg(34)),
    "G35" => (kind = :complex, K = 1, group = () -> Chevie.crg(35)),
    "G36" => (kind = :complex, K = 1, group = () -> Chevie.crg(36)),
    "G37" => (kind = :complex, K = 2, group = () -> Chevie.crg(37)),
    "dual_H4" => (kind = :dual, K = 1, group = () -> Chevie.coxgroup(:H, 4)),
    "dual_D8" => (kind = :dual, K = 1, group = () -> Chevie.coxgroup(:D, 8)),
    "classic_E8" => (kind = :classic, K = 0, group = () -> Chevie.coxgroup(:E, 8)),
)

ngens(W) = length(Chevie.gens(W))

monoid(kind, W) = kind === :classic ? Chevie.BraidMonoid(W) :
                  kind === :dual ? Chevie.DualBraidMonoid(W; c = collect(1:ngens(W))) :
                  Chevie.DualBraidMonoid(W)

germ_p2(kind, W) = kind === :classic ? GarsideHomology.CoxeterGerm(Chevie.cartan(W)) :
                   kind === :dual ? dual_germ(W; c = collect(1:ngens(W))) : dual_germ(W)

resolution(G) = lcm_resolution(LCMProjectionData(CellTree(G)))

# Warm-up: the same functions on small germs of the same types.
function warm_up(kind, path, K)
    if kind === :classic
        path == "P1" && germ(Chevie.BraidMonoid(Chevie.coxgroup(:A, 3)))
        resolution(GarsideHomology.CoxeterGerm(Chevie.cartan(Chevie.coxgroup(:A, 8))))   # rank 8, lazy
        return nothing
    end
    W = kind === :complex ? Chevie.crg(4) : Chevie.coxgroup(:I, 2, 5)
    resolution(path == "P1" ? germ(monoid(kind, W)) : germ_p2(kind, W))
    if K == 2                                    # descents on two words: 65 atoms
        W2 = Chevie.coxgroup(:I, 2, 65)
        resolution(path == "P1" ? germ(monoid(:dual, W2)) : germ_p2(:dual, W2))
    end
    return nothing
end

gio(x) = round(x / 2^30; digits = 3)

# A step: @timed, then the peak and the live memory after a full GC.
function measure_step(f)
    st = @timed f()
    rss = Sys.maxrss()
    GC.gc(true)
    live = Base.gc_live_bytes()
    return (value = st.value, time = st.time, gc = st.gctime, compile = st.compile_time,
            bytes = st.bytes, rss = rss, live = live)
end

const STEPS = ["groupe", "monoide", "germe", "arbre", "diff"]

function append_csv(path, header, row)
    new = !isfile(path)
    open(path, "a") do io
        new && println(io, join(header, ','))
        println(io, join((replace(string(x), ',' => ';') for x in row), ','))
    end
end

git(args...) = try readchomp(Cmd(["git", "-C", dirname(@__DIR__), args...])) catch; "?" end

function main()
    length(ARGS) >= 4 || error("Usage: julia --project=benchmark benchmark/construction.jl <machine> <case> <P1|P2> <round>")
    machine, id, path, round_ = ARGS[1], ARGS[2], ARGS[3], parse(Int, ARGS[4])
    haskey(CASES, id) || error("unknown case $id; cases: $(join(sort(collect(keys(CASES))), ", "))")
    path in ("P1", "P2") || error("path $path: P1 or P2")
    case = CASES[id]
    warm_up(case.kind, path, case.K)
    GC.gc(true)
    rss0 = Sys.maxrss()
    s = Dict{String, Any}()
    s["groupe"] = measure_step(case.group)
    W = s["groupe"].value
    if path == "P1"
        s["monoide"] = measure_step(() -> monoid(case.kind, W))
        M = s["monoide"].value
        s["germe"] = measure_step(() -> germ(M))
    else
        s["germe"] = measure_step(() -> germ_p2(case.kind, W))
    end
    G = s["germe"].value
    s["arbre"] = measure_step(() -> CellTree(G))
    T = s["arbre"].value
    s["diff"] = measure_step(function ()
        P = LCMProjectionData(T)
        C = lcm_resolution(P)
        return (C, GarsideHomology._nchains(P.chains))
    end)
    C, chains = s["diff"].value
    terms = sum(GarsideHomology.Complexes.SparseArrays.nnz(differential(C, n)) for n in 0:(last(range(C)) - 1); init = 0)
    fp = applicable(simples, G) ? sprint(show, GarsideHomology.fingerprint(G)) : "germe paresseux"
    header = ["machine", "julia", "threads", "commit_gauss", "modifs_gauss", "cas", "chemin", "tour",
              "cellules", "termes", "chaines", "maxrss_debut_Gio"]
    row = Any[machine, string(VERSION), Threads.nthreads(), git("rev-parse", "--short", "HEAD"),
              isempty(git("status", "--porcelain", "--untracked-files=no")) ? "non" : "oui",
              id, path, round_, ncells(T), terms, chains, gio(rss0)]
    for k in STEPS
        append!(header, ["t_$(k)_s", "gc_$(k)_s", "compil_$(k)_s", "alloc_$(k)_Gio", "maxrss_$(k)_Gio", "vivant_$(k)_Gio"])
        if haskey(s, k)
            x = s[k]
            append!(row, [round(x.time; sigdigits = 5), round(x.gc; sigdigits = 4), round(x.compile; sigdigits = 3),
                          gio(x.bytes), gio(x.rss), gio(x.live)])
        else
            append!(row, fill("", 6))
        end
    end
    append!(header, ["empreinte_germe", "date"])
    append!(row, [fp, Libc.strftime("%Y-%m-%dT%H:%M:%S", time())])
    mkpath(OUT_DIR)
    append_csv(joinpath(OUT_DIR, "construction_$(machine).csv"), header, row)
    println("$id $path, tour $round_ : ",
            path == "P1" ? "monoïde $(round(s["monoide"].time; digits = 1)) s, " : "",
            "germe $(round(s["germe"].time; digits = 2)) s, arbre $(round(s["arbre"].time; digits = 2)) s, ",
            "différentielle $(round(s["diff"].time; digits = 2)) s, pic $(gio(s["diff"].rss)) Gio, ",
            "$(ncells(T)) cellules, $terms termes")
end

main()
