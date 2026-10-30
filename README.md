# Romanoff's constant: C₄₅ < 0.489481743486002

A certificate for an upper bound on [Romanoff's constant](https://teorth.github.io/optimizationproblems/constants/45a.html)
(C₄₅ in Tao's optimization-constants repository): the upper density of

    A = {n odd : n = p + 2^k, p prime, k ≥ 0}

satisfies **d̄(A) < 0.489481743486002**. The certified rational sum is
S = 0.48948174348600104752999373053459010942525432665839… (numerator and denominator in
`certs/c45-54/summary.json`).

| bound | source |
|---|---|
| 0.490249407811155 | Griego, 36 primes |
| 0.490180063290061 | Yoo, 38 primes |
| 0.489596115074071 | Patankar, open PR [teorth/optimizationproblems#281](https://github.com/teorth/optimizationproblems/pull/281), 21-prime base |
| 0.489520050199154 | this repository, 51 primes (`cert-51.json`) |
| 0.489507858991149 | this repository, 52 primes (`cert-52.json`) |
| **0.489481743486002** | this repository, 54 primes (`cert.json`) |

## Verify

```bash
./verify.sh --stored        # check the stored histograms (about 20 s)
./verify.sh cert-51.json    # rebuild the 51-prime histogram from the residue conditions, then check (~17 min)
./verify.sh                 # rebuild the 54-prime histogram, then check (~1 h with THREADS=2, ~4 GB)
```

A C++20 compiler and Python 3.10+ (standard library only). `EXPECTED_OUTPUT*.txt` are the verifier's outputs.
The rebuilt histograms must hash to (sha256 of the uncompressed text):

- 54 primes: 1,434,945 bins + 1 aggregate line, `e011429b8f49cd362db1f8a3273fe5ceeff97d1067af1e97bfbbdf2bddfaf91c`
- 52 primes: 592,715 entries, `3c89e3e0a5409e63ab73766a1dbdfbd5aca143a6d31145465797af2632f0dff9`
- 51 primes: 56,798 entries, `67943dd51e6df1e671741f8ad72b22283899a8c976d93cd9549472a003857b6f`

## The inequality

For a finite set P of odd primes, M = ∏P, T = lcm ord_q(2) (Erdős' obstruction argument, as used by
Habsieger–Roblot, Chen–Dai–Li, Griego, Yoo):

    d̄(A) ≤ Σ_{a mod M} min( 1/(2M), ν(a) / (T φ(M) log 2) ),   ν(a) = #{r mod T : gcd(a − 2^r, M) = 1}.

P is split into a *base* (exact histogram of ν₀ over a mod M₀) and *clusters* (primes of one odd order,
coprime to the base period and to each other, so ν = ν₀ · ∏ h_d by CRT). log 2 is replaced by the rational
lower bound 2 Σ_{j<40} 1/((2j+1) 3^(2j+1)), which only increases the sum.

## The set (54 primes)

Base, 29 primes, orders dividing T₀ = 12252240 = 2⁴·3²·5·7·11·13·17:

- 3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73, 257 (Yoo's base)
- 43, 127, 241, 113, 109, 97, 151 (added in #281)
- **8191, 53, 157** (orders 13, 52, 52), **23, 89, 397, 683** (orders 11, 11, 44, 22) and **137** (order 68), new here

Clusters, 25 primes: 23: 47, 178481 | 29: 233, 1103, 2089 | 37: 223 | 43: 431, 9719 | 47: 2351, 4513 |
53: 6361 | 73: 439 | 79: 2687 | 83: 167 | 113: 3391 | 131: 263 | 179: 359, 1433 | 191: 383 | 233: 1399 |
239: 479, 1913, 5737, 176383 | 251: 503 (Griego's and Yoo's clusters, without orders 11 and 13).

Why the new primes were missed: 8191 (order 13) and 23, 89 (order 11) sat in clusters, which keeps 13 and 11
out of the base period. Moving them into the base leaves the sum unchanged (same CRT factorisation) and makes
primes of order 2p and 4p available. Measured gains: {53, 157} 7.61e-5, 397 1.22e-5, {683, 137} 2.61e-5.

The 51- and 52-prime certificates are subsets: without 23, 89, 397, 683, 137 in the base, and without 683, 137.

## How the base histogram is computed

Grouping the base primes by the part of their order outside 12 (2-part beyond 4, 3-part beyond 3, each prime
≥ 5) splits r mod T₀ into r mod 12 and one private coordinate per group, and a prime only constrains r mod 12
and its group's coordinate. Hence

    ν₀ = Σ_{c mod 12} ∏_G α_G(c),   α_G(c) = #{x mod 12·P_G : x ≡ c (12), x admissible for G},

with independent α_G. Each group's law is an exact DP; partial products are merged up to the affine maps of
Z/12. Groups whose orders are prime to 3 see only r mod 4: with σ_s = Σ_{c ≡ s (4)} of the product of the other
groups, ν₀ = Σ_{s mod 4} σ_s β_s, where σ ranges over its orbits under the affine maps of Z/4 and β over every
element of the law of the other groups. Counts are 256-bit (M₀ ≈ 10⁵²).

For the 54-prime certificate the pairs (σ, β) are summed per weight of β, and a pair is not enumerated when
its ν₀ is certainly above C = 1,538,836 (the entry sum of σ times the least entry of β already exceeds C).
Those classes go to one aggregate line `# above C count Σν₀ Σν₀²`, computed exactly from suffix sums. For
ν₀ > C the term is 1/(2M) for every cluster multiplier, which the verifier checks ((C+1)·m_min ≥ thr), so the
aggregate changes nothing in S. On the 51- and 52-prime sets this path gives histograms byte-identical to the
plain one without the cut, and with the cut the same bins up to C and the same count, Σν₀, Σν₀² above it.
`dp_seed2.cpp --self-test` checks every mode against direct enumeration over a mod M₀ on small prime sets.

`verify_cert.py` does not trust the generator: it checks the histogram's total (= M₀) and its first two moments,
aggregate included, against values computed directly from the residue conditions, checks its own pipeline
against the bound computed from its definition on a small prime set (and that a cut below the trivial regime is
rejected), and evaluates S with `Fraction`s.

## Cross-checks

- Griego's and Yoo's seed histograms are reproduced byte for byte, and their sums to all 100 printed digits.
- The 21-prime base histogram of #281 is reproduced byte for byte (sha256 `238174d0…fcc41d`) by two different
  modes of the generator; this is an independent replay of that result.
- The 51-prime histogram was rebuilt from this package in a clean copy and matches.

## Credits

The inequality and the base/cluster certificate format are Sebastian Griego's
([c45-romanoff-certificate](https://github.com/sebastian-griego/c45-romanoff-certificate)); Andrew Yoo improved
the prime set ([romanoff v3](https://github.com/andrew-yoo/romanoff/tree/v3)); V. Patankar's open PR
[#281](https://github.com/teorth/optimizationproblems/pull/281) extended the base to 21 primes, which this base
contains. The generator and verifier here were written independently of theirs.

Found and prepared with Claude (Anthropic) agents under Rustam Salavatov's direction (search, generator,
verifier and text). The value is recomputed by the code in this repository; it has not been verified line by
line by a human expert.

MIT license.
