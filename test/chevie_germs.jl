# Comparison of a germ built from Chevie with Chevie itself (J1.c, c3): a tool of the
# tests of GarsideHomology, also included by the scripts of CatLab (j1/c3_chevie.jl).

using GarsideHomology
import Chevie

"""
    chevie_elements(G, M)

The element of Chevie of each simple of `G`, built from `M`: the product of the atoms of
`M` along `word(G, s)`, at the position `s + 1`.
"""
chevie_elements(G, M) =
    [foldl(*, (M.atoms[a] for a in GarsideHomology.word(G, s)); init = one(M)) for s in GarsideHomology.simples(G)]

# The pairs (s, a) whose lcm is compared: all of them if there are at most `max_pairs`,
# otherwise a reproducible sample of `sample` pairs and all the pairs where s is an atom
# or Δ.
function _lcm_pairs(n, N; max_pairs, sample, seed)
    n * N <= max_pairs && return [(Int32(s), a) for s in 0:(n - 1) for a in 1:N], true
    chosen = Set{Tuple{Int32, Int}}()
    for s in vcat(1:N, n - 1), a in 1:N
        push!(chosen, (Int32(s), a))
    end
    x = UInt64(seed)
    for _ in 1:sample
        x = x * 0x5851f42d4c957f2d + 0x14057b7ef767814f
        push!(chosen, (Int32((x >> 16) % n), Int((x >> 48) % N) + 1))
    end
    return sort!(collect(chosen)), false
end

"""
    compare_with_chevie(G, M; max_pairs = 20_000, sample = 2_000, seed = 1)

Compare the germ `G = TableGerm(M)` with the monoid `M` of Chevie: the elements of the
simples are distinct and as many as in Chevie (the simples of `left_divisors` for a dual
monoid, |W| for a braid monoid); for each simple, its word (`Chevie.word`, the smallest
reduced word), its length and its left descents (`Chevie.leftdescents`); for pairs (s, a)
chosen by `_lcm_pairs`, the lcm with the atom and its complement (`Chevie.rightlcmc`).
Returns `(failures, pairs, complete)`; the first failures are reported by warnings.
"""
function compare_with_chevie(G, M; max_pairs = 20_000, sample = 2_000, seed = 1)
    failures = 0
    fail(msg) = (failures += 1; failures <= 20 && @warn msg; nothing)
    S = GarsideHomology.simples(G)
    n, N = length(S), natoms(G)
    N == length(M.atoms) || fail("$N atoms, Chevie has $(length(M.atoms))")
    elt = chevie_elements(G, M)
    length(unique(elt)) == n || fail("$(n - length(unique(elt))) simples share their element of Chevie")
    expected = M isa Chevie.DualBraidMonoid ? sum(length, Chevie.left_divisors(M, M.δ)) :
                                             prod(Chevie.degrees(M.W))
    n == expected || fail("$n simples, Chevie has $expected")
    for s in S
        e = elt[s + 1]
        w = Chevie.word(M, e)
        w == GarsideHomology.word(G, s) || fail("word of $s: $(GarsideHomology.word(G, s)), Chevie: $w")
        length(G, s) == length(w) || fail("length of $s: $(length(G, s)), Chevie: $(length(w))")
        Chevie.leftdescents(M, e) == collect(GarsideHomology.left_descents(G, s)) ||
            fail("left descents of $s: $(collect(GarsideHomology.left_descents(G, s))), " *
                 "Chevie: $(Chevie.leftdescents(M, e))")
    end
    pairs, complete = _lcm_pairs(n, N; max_pairs, sample, seed)
    for (s, a) in pairs
        m, (τ, _) = Chevie.rightlcmc(M, elt[s + 1], M.atoms[a])
        r = GarsideHomology.lcm_atom(G, s, a)
        elt[r.lcm + 1] == m && elt[r.complement + 1] == τ ||
            fail("lcm of $s and atom $a: $r, which is not the lcm of Chevie")
    end
    return (failures = failures, pairs = length(pairs), complete = complete)
end
