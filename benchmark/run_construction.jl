# benchmark/run_construction.jl — campaign of timings of the construction of the lcm
# resolution (J1.d, task d7), protocol of the task C1 of J0
#
# Usage, from the root of GarsideHomology.jl, the machine free (nothing else running):
#     julia --project=benchmark -e "using Pkg; Pkg.instantiate()"      # once
#     julia --project=benchmark benchmark/run_construction.jl N9           # 7 cases × P1, P2 × 3 rounds
#     julia --project=benchmark benchmark/run_construction.jl N9 G37 G34   # chosen cases
#     julia --project=benchmark benchmark/run_construction.jl N9 j0=<path of temps_N9.csv of J0>
# Runs construction.jl in a fresh Julia session for each round (--threads=auto, written in
# the CSV), in random order, and skips the rounds already in benchmark/out/construction_<machine>.csv:
# the campaign can be interrupted and run again. Log of each round in benchmark/out/logs/.
# At the end, the minima by case and path, beside those of J0 (j0/out/temps_<machine>.csv
# of CatLab, by default ../CatLab/j0/out/ beside GarsideHomology.jl).

# Only Base: no Printf, no Random (Pkg did not take them as standard libraries on N9).

const MACHINE = isempty(ARGS) ? error("Usage: julia --project=benchmark benchmark/run_construction.jl <machine> [cases] [j0=<csv>]") : ARGS[1]
const REFERENCE = ["G35", "G36", "dual_H4", "G37", "dual_D8", "classic_E8", "G34"]
const OPTIONS = filter(a -> startswith(a, "j0="), ARGS[2:end])
const CHOSEN = filter(a -> !startswith(a, "j0="), ARGS[2:end])
const CASES = isempty(CHOSEN) ? REFERENCE : CHOSEN
const PATHS = ["P1", "P2"]
const ROUNDS = 3
const OUT_DIR = joinpath(@__DIR__, "out")
const LOG_DIR = joinpath(OUT_DIR, "logs")
const CSV = joinpath(OUT_DIR, "construction_$(MACHINE).csv")
const J0_CSV = isempty(OPTIONS) ? normpath(joinpath(@__DIR__, "..", "..", "CatLab", "j0", "out", "temps_$(MACHINE).csv")) :
               String(chop(OPTIONS[end]; head = 3, tail = 0))

function read_csv(path)
    isfile(path) || return String[], Vector{Vector{String}}()
    lines = readlines(path)
    return split(lines[1], ','), [split(l, ',') for l in lines[2:end] if !isempty(l)]
end

column(header, name) = findfirst(==(name), header)

function done()
    h, rows = read_csv(CSV)
    isempty(h) && return Set{Tuple{String, String, Int}}()
    ic, ip, it = column(h, "cas"), column(h, "chemin"), column(h, "tour")
    return Set((String(r[ic]), String(r[ip]), parse(Int, r[it])) for r in rows)
end

# Fisher–Yates, with rand of Base.
function shuffle_rounds!(v)
    for i in length(v):-1:2
        j = rand(1:i)
        v[i], v[j] = v[j], v[i]
    end
    return v
end

function campaign()
    mkpath(LOG_DIR)
    finished = done()
    todo = shuffle_rounds!([(c, p, r) for c in CASES, p in PATHS, r in 1:ROUNDS if !((c, p, r) in finished)])
    println(length(todo), " tours à faire ; CSV : ", CSV)
    script = joinpath(@__DIR__, "construction.jl")
    for (k, (c, p, r)) in enumerate(todo)
        log = joinpath(LOG_DIR, "construction_$(c)_$(p)_$(r).log")
        print("[$k/$(length(todo))] $c $p tour $r … ")
        flush(stdout)
        cmd = `$(Base.julia_cmd()) --project=$(@__DIR__) --threads=auto $script $MACHINE $c $p $r`
        t = @elapsed ok = open(io -> success(pipeline(cmd; stdout = io, stderr = io)), log, "w")
        println(ok ? "fait" : "ÉCHEC (voir $log)", " en ", round(t; digits = 1), " s")
    end
end

minimum_or_nan(xs) = isempty(xs) ? NaN : minimum(xs)

function summary()
    h, rows = read_csv(CSV)
    isempty(h) && return
    get_(r, name) = (i = column(h, name); i === nothing || isempty(r[i]) ? NaN : parse(Float64, r[i]))
    hj, rj = read_csv(J0_CSV)
    println()
    isempty(hj) ? println("Temps de J0 absents ($J0_CSV).") : println("Temps de J0 : ", J0_CSV)
    println("Minima sur les tours (s, Gio) ; J0 : arbre + différentielle de gauss.jl, Gauss_a")
    widths = [11, 9, 9, 7, 8, 8, 9, 8, 7, 8, 9, 8, 7, 9, 9]
    titles = ["cas", "J0 arbre", "J0 diff", "J0 pic", "P1 mon.", "germe", "arbre", "diff", "pic",
              "P2 germe", "arbre", "diff", "pic", "gain", "gain+g"]
    println("gain : (arbre + différentielle de J0) / (arbre + différentielle de P1) ; ",
            "gain+g : la même chose avec le germe de P1 au dénominateur")
    println(rpad(titles[1], widths[1]), " ", join((lpad(t, w) for (t, w) in zip(titles[2:end], widths[2:end])), " "))
    for c in CASES
        sel(p) = [r for r in rows if r[column(h, "cas")] == c && r[column(h, "chemin")] == p]
        m(p, name) = minimum_or_nan([get_(r, name) for r in sel(p)])
        j0_rows = isempty(hj) ? [] : [r for r in rj if r[column(hj, "cas")] == c]
        j(name) = minimum_or_nan([parse(Float64, r[column(hj, name)]) for r in j0_rows])
        j0_time = j("t_arbre_s") + j("t_diff_s")
        gain = j0_time / (m("P1", "t_arbre_s") + m("P1", "t_diff_s"))
        gain_germ = j0_time / (m("P1", "t_germe_s") + m("P1", "t_arbre_s") + m("P1", "t_diff_s"))
        values = [j("t_arbre_s"), j("t_diff_s"), j("maxrss_fin_Gio"),
                  m("P1", "t_monoide_s"), m("P1", "t_germe_s"), m("P1", "t_arbre_s"), m("P1", "t_diff_s"), m("P1", "maxrss_diff_Gio"),
                  m("P2", "t_germe_s"), m("P2", "t_arbre_s"), m("P2", "t_diff_s"), m("P2", "maxrss_diff_Gio"),
                  gain, gain_germ]
        # Times and gains with three significant digits, peaks with two decimals.
        peak = [false, false, true, false, false, false, false, true, false, false, false, true, false, false]
        show_value(x, p) = isnan(x) ? "-" : p ? string(round(x; digits = 2)) : string(round(x; sigdigits = 3))
        println(rpad(c, widths[1]), " ",
                join((lpad(show_value(x, p), w) for (x, p, w) in zip(values, peak, widths[2:end])), " "))
    end
end

campaign()
summary()
