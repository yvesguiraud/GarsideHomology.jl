# Verification of a germ whose simples are all numbered, and counts of multichains.

# Calls f(ν) for each multiple ν of ω (ω included), by layers of increasing length along
# the right covers; `marks` (one entry per simple) and `tag` avoid visiting twice.
function _foreach_multiple(f, G::AbstractGerm, ω::Int32, marks::Vector{Int32}, tag::Int32,
                           layer::Vector{Int32}, next::Vector{Int32})
    N = natoms(G)
    empty!(layer)
    push!(layer, ω)
    marks[ω + 1] = tag
    f(ω)
    while !isempty(layer)
        empty!(next)
        for m in layer, a in 1:N
            p = product(G, m, a)
            (p === nothing || marks[p + 1] == tag) && continue
            marks[p + 1] = tag
            push!(next, p)
            f(p)
        end
        layer, next = next, layer
    end
    return nothing
end

# Calls f(σ, ν) for each pair σ ≼ ν of simples.
function _foreach_interval(f, G::AbstractGerm)
    marks = zeros(Int32, nsimples(G))
    layer, next = Int32[], Int32[]
    for σ in simples(G)
        _foreach_multiple(ν -> f(σ, ν), G, σ, marks, σ + Int32(1), layer, next)
    end
    return nothing
end

"""
    count_multichains(G::AbstractGerm, m) -> Int

The number of multichains σ₁ ≼ σ₂ ≼ ⋯ ≼ σₘ of `m` simples of `G` (repetitions allowed),
whose simples must all be numbered: the number of simples for `m = 1`, of pairs σ ≼ ν for
`m = 2`. For the dual braid monoid of a well-generated group, this is the Fuss–Catalan
number Π (dᵢ + m·h)/dᵢ. Computed by m − 1 passes over the intervals, in O(n) memory; an
`OverflowError` is thrown rather than a wrong count.
"""
function count_multichains(G::AbstractGerm, m::Integer)
    m >= 1 || throw(ArgumentError("count_multichains: m = $m < 1"))
    f = ones(Int, nsimples(G))
    for _ in 2:m
        f = _multichain_step(G, f)
    end
    return foldl(Base.checked_add, f; init = 0)
end

function _multichain_step(G::AbstractGerm, f::Vector{Int})
    g = zeros(Int, length(f))
    _foreach_interval(G) do σ, ν
        g[ν + 1] = Base.checked_add(g[ν + 1], f[σ + 1])
    end
    return g
end

"""
    verify(G::AbstractGerm; full = false)

Check the consistency of the tables of a germ whose simples are all numbered (a tabulated
germ, a germ of the tests, a germ read from a file), and throw an `ArgumentError` naming
the first defect. Checked, in O(n·N) operations of the interface:
- the identity and the atoms (lengths 0 and 1, words, descents);
- the left descents of each simple are exactly the atoms that divide it, with quotients
  of length one less;
- each simple product by an atom adds 1 to the length and divides back;
- for each simple ω and atom a, `lcm_atom` gives a common multiple ω·τ = ω ∨ a, equal to
  (ω, 1) when a ≼ ω;
- a single simple Δ without right cover, which every simple divides;
- the canonical order (by length, then smallest word) and the lexicographic comparison.
With `full = true`, also the lattice property, in O(#intervals · N): each ω ∨ a divides
every common multiple of ω and a.
"""
function verify(G::AbstractGerm; full::Bool = false)
    fail(msg) = throw(ArgumentError("verify: " * msg))
    S = simples(G)
    n, N = length(S), natoms(G)
    S == 0:(n - 1) || fail("the simples are not 0:$(n - 1)")
    length(G, 0) == 0 && isempty(left_descents(G, 0)) || fail("the identity has a length or a descent")
    for a in 1:N
        length(G, a) == 1 && word(G, a) == [a] && collect(left_descents(G, a)) == [a] ||
            fail("atom $a: length, word or descents")
    end
    maximal = [s for s in S if all(a -> product(G, s, a) === nothing, 1:N)]
    length(maximal) == 1 || fail("$(length(maximal)) simples without right cover")
    Δ = only(maximal)
    previous = (0, Int32[])
    for s in S
        ℓ = length(G, s)
        D = left_descents(G, s)
        for a in 1:N
            q = ldiv(G, a, s)
            (q !== nothing) == (a in D) || fail("left_descents($s) and the division by atom $a")
            q === nothing || (length(G, q) == ℓ - 1 && product(G, a, q) == s) ||
                fail("quotient of $s by atom $a")
            p = product(G, s, a)
            p === nothing || (length(G, p) == ℓ + 1 && ldiv(G, s, p) == a) ||
                fail("product of $s by atom $a")
            r = lcm_atom(G, s, a)
            product(G, s, r.complement) == r.lcm && is_left_divisor(G, a, r.lcm) ||
                fail("lcm_atom($s, $a) = $r is not a common multiple")
            a in D && !(r.lcm == s && iszero(r.complement)) && fail("lcm_atom($s, $a) with $a ≼ $s")
        end
        is_left_divisor(G, s, Δ) || fail("$s does not divide Δ = $Δ")
        w = word(G, s)
        length(w) == ℓ || fail("word($s) of length $(length(w)) ≠ $ℓ")
        if s > 0
            previous < (ℓ, w) || fail("canonical order at $s")
            isless_lex(G, s - 1, s) == isless(previous[2], w) || fail("isless_lex($(s - 1), $s)")
        end
        previous = (ℓ, w)
    end
    full && _verify_lattice(G)
    return nothing
end

# Each ω ∨ a divides every multiple of ω that a divides.
function _verify_lattice(G::AbstractGerm)
    marks = zeros(Int32, nsimples(G))
    layer, next = Int32[], Int32[]
    for ω in simples(G)
        _foreach_multiple(G, ω, marks, ω + Int32(1), layer, next) do m
            for a in left_descents(G, m)
                is_left_divisor(G, lcm_atom(G, ω, a).lcm, m) ||
                    throw(ArgumentError("verify: lcm_atom($ω, $a) does not divide the common multiple $m"))
            end
        end
    end
    return nothing
end
