# Romanoff's constant: C₄₅ < 0.489507858991149

A certificate for an upper bound on [Romanoff's constant](https://teorth.github.io/optimizationproblems/constants/45a.html)
(C₄₅ in Tao's optimization-constants repository): the upper density of

    A = {n odd : n = p + 2^k, p prime, k ≥ 0}

satisfies **d̄(A) < 0.489507858991149**. The certified rational sum is
S = 0.48950785899114820764182239975315014324271757802496723921615… (numerator and denominator in
`certs/c45-52/summary.json`).

| bound | source |
|---|---|
| 0.490249407811155 | Griego, 36 primes |
| 0.490180063290061 | Yoo, 38 primes |
| 0.489596115074071 | Patankar, open PR [teorth/optimizationproblems#281](https://github.com/teorth/optimizationproblems/pull/281), 21-prime base |
| 0.489520050199154 | this repository, 51 primes (`cert-51.json`) |
| **0.489507858991149** | this repository, 52 primes (`cert.json`) |

## Verify

```bash
./verify.sh --stored        # check the stored histograms (seconds)
./verify.sh cert-51.json    # rebuild the 51-prime histogram from the residue conditions, then check (~17 min)
./verify.sh                 # rebuild the 52-prime histogram, then check (~1.5 h with THREADS=2, ~3.6 GB)
```

A C++20 compiler and Python 3.10+ (standard library only). `EXPECTED_OUTPUT.txt` and
`EXPECTED_OUTPUT-51.txt` are the verifier's outputs. The rebuilt histograms must hash to

- 52 primes: 592,715 entries, sha256 `3c89e3e0a5409e63ab73766a1dbdfbd5aca143a6d31145465797af2632f0dff9`
- 51 primes: 56,798 entries, sha256 `67943dd51e6df1e671741f8ad72b22283899a8c976d93cd9549472a003857b6f`

## The inequality

For a finite set P of odd primes, M = ∏P, T = lcm ord_q(2) (Erdős' obstruction argument, as used by
Habsieger–Roblot, Chen–Dai–Li, Griego, Yoo):

    d̄(A) ≤ Σ_{a mod M} min( 1/(2M), ν(a) / (T φ(M) log 2) ),   ν(a) = #{r mod T : gcd(a − 2^r, M) = 1}.

P is split into a *base* (exact histogram of ν₀ over a mod M₀) and *clusters* (primes of one odd order,
coprime to the base period and to each other, so ν = ν₀ · ∏ h_d by CRT). log 2 is replaced by the rational
lower bound 2 Σ_{j<40} 1/((2j+1) 3^(2j+1)), which only increases the sum.

## The set (52 primes)

Base, 27 primes, orders dividing T₀ = 720720 = 2⁴·3²·5·7·11·13:

- 3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73, 257 (Yoo's base)
- 43, 127, 241, 113, 109, 97, 151 (added in #281)
- **8191, 53, 157** (orders 13, 52, 52) and **23, 89, 397** (orders 11, 11, 44) (new here)

Clusters, 25 primes: 23: 47, 178481 | 29: 233, 1103, 2089 | 37: 223 | 43: 431, 9719 | 47: 2351, 4513 |
53: 6361 | 73: 439 | 79: 2687 | 83: 167 | 113: 3391 | 131: 263 | 179: 359, 1433 | 191: 383 | 233: 1399 |
239: 479, 1913, 5737, 176383 | 251: 503 (Griego's and Yoo's clusters, without orders 11 and 13).

Why the new primes were missed: 8191 (order 13) and 23, 89 (order 11) sat in clusters, which keeps 13 and 11
out of the base period. Moving them into the base leaves the sum unchanged (same CRT factorisation) and makes
primes of order 4·13 and 4·11 available. {53, 157} lower the bound by 7.61e-5 and 397 by 1.22e-5.

The 51-prime certificate is the same set without 23, 89, 397 in the base (23, 89 stay a cluster).

## How the base histogram is computed

Grouping the base primes by the part of their order outside 12 (2-part beyond 4, 3-part beyond 3, each prime
≥ 5) splits r mod T₀ into r mod 12 and one private coordinate per group, and a prime only constrains r mod 12
and its group's coordinate. Hence

    ν₀ = Σ_{c mod 12} ∏_G α_G(c),   α_G(c) = #{x mod 12·P_G : x ≡ c (12), x admissible for G},

with independent α_G. Each group's law is an exact DP; partial products are merged up to the affine maps of
Z/12 (the joint law is invariant under them). Groups whose orders are prime to 3 see only r mod 4: with
σ_s = Σ_{c ≡ s (4)} of the product of the other groups, ν₀ = Σ_{s mod 4} σ_s β_s, and σ orbits under the
affine maps of Z/4 are paired with every element of the law of β. Final counts are 256-bit (M₀ ≈ 10⁴⁷).
`dp_seed2.cpp --self-test` checks every mode against direct enumeration over a mod M₀ on small prime sets.

`verify_cert.py` does not trust the generator: it checks the histogram's total (= M₀) and its first two moments
against values computed directly from the residue conditions, checks its own pipeline against the bound
computed from its definition on a small prime set, and evaluates S with `Fraction`s.

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
