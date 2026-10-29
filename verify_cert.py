"""Exact verifier for a Romanoff upper-density certificate (base histogram + clusters).

    python3 verify_cert.py CERT.json SEED_HISTOGRAM.txt [--self-test] [--json-out OUT]

CERT.json: {"base_order": [...], "clusters": {"d": [q, ...], ...}}. The seed histogram is
the output of dp_seed for base_order. Only the standard library; every step of the bound is
integer or Fraction arithmetic, decimals are for display.

Bound (Erdos' obstruction argument, as in Habsieger-Roblot, Chen-Dai-Li, Griego, Yoo):
for P = base U clusters, M = prod P, T = lcm ord_q(2),
    d(A) <= sum_{a mod M} min(1/(2M), nu(a)/(T phi(M) log 2)),
    nu(a) = #{r mod T : gcd(a - 2^r, M) = 1}.
Cluster orders are pairwise coprime and coprime to T0 = lcm of base orders, so by CRT
nu(a) = nu0(a mod M0) * prod_d h_d(a mod Q_d), each h_d from cluster_poly.
log 2 is replaced by a rational lower bound, which only increases the sum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from bisect import bisect_right
from collections import Counter
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import product
from math import gcd, lcm, prod
from pathlib import Path


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def ord2(q: int) -> int:
    x, d = 2 % q, 1
    while x != 1:
        x, d = x * 2 % q, d + 1
    return d


def log2_lower(terms: int = 40) -> Fraction:
    # log 2 = 2 sum_j 1/((2j+1) 3^(2j+1)); partial sums are lower bounds
    return sum(Fraction(2, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(terms))


def cluster_poly(d: int, primes: list[int]) -> dict[int, int]:
    dist = {0: 1}
    for q in primes:
        require(ord2(q) == d, f"ord_{q}(2) != {d}")
        new: Counter[int] = Counter()
        for s, c in dist.items():
            new[s] += c * (q - d)
            if s:
                new[s] += c * s
            if s < d:
                new[s + 1] += c * (d - s)
        dist = dict(new)
    poly = {d - s: c for s, c in dist.items()}
    require(sum(poly.values()) == prod(primes), "cluster weight mismatch")
    return poly


def cluster_poly_brute(d: int, primes: list[int]) -> dict[int, int]:
    pw = [[pow(2, r, q) for r in range(d)] for q in primes]
    out: Counter[int] = Counter()
    for a in product(*(range(q) for q in primes)):
        out[sum(all(a[i] != pw[i][r] for i in range(len(primes))) for r in range(d))] += 1
    return dict(out)


def atoms(clusters: dict[int, list[int]]) -> list[tuple[int, int]]:
    mult: Counter[int] = Counter({1: 1})
    for d, qs in clusters.items():
        poly = cluster_poly(d, qs)
        new: Counter[int] = Counter()
        for m, w in mult.items():
            for h, c in poly.items():
                new[m * h] += w * c
        mult = new
    return sorted(mult.items())


def upper_sum(hist: dict[int, int], base: list[int], clusters: dict[int, list[int]],
              split: int | None = None, terms: int = 40) -> Fraction:
    T0 = 1
    for q in base:
        T0 = lcm(T0, ord2(q))
    allp = base + [q for qs in clusters.values() for q in qs]
    require(len(set(allp)) == len(allp), "duplicate prime")
    require(all(is_prime(q) and q > 2 for q in allp), "not an odd prime")
    ds = list(clusters)
    for i, d in enumerate(ds):
        require(gcd(d, T0) == 1, f"cluster order {d} shares a factor with T0={T0}")
        for e in ds[i + 1:]:
            require(gcd(d, e) == 1, f"cluster orders {d}, {e} not coprime")
    require(sum(hist.values()) == prod(base), "seed histogram total != prod(base)")
    require(all(0 <= nu <= T0 for nu in hist), "nu out of range")

    M = prod(allp)
    phi = prod(q - 1 for q in allp)
    T = 1
    for q in allp:
        T = lcm(T, ord2(q))

    ds = sorted(clusters, key=lambda d: len(clusters[d]))
    if split is None:
        best = None
        for s in range(len(ds) + 1):
            a = prod(len(clusters[d]) + 1 for d in ds[:s])
            b = prod(len(clusters[d]) + 1 for d in ds[s:])
            cost = a * len(hist) * 20 + b
            if best is None or cost < best:
                best, split = cost, s
    A1 = atoms({d: clusters[d] for d in ds[:split]})
    A2 = atoms({d: clusters[d] for d in ds[split:]})
    f2 = [m for m, _ in A2]
    pw, pwm = [0], [0]
    for m, w in A2:
        pw.append(pw[-1] + w)
        pwm.append(pwm[-1] + w * m)
    require(pw[-1] * sum(w for _, w in A1) == prod(q for qs in clusters.values() for q in qs),
            "multiplier weight mismatch")

    l2 = log2_lower(terms)
    # nu/(T phi l2) <= 1/(2M)  <=>  nu <= thr
    thr = Fraction(T * phi, 2 * M) * l2
    tn, td = thr.numerator, thr.denominator
    small = large = 0
    for nu0, delta in hist.items():
        if nu0 == 0:
            continue
        for m1, w1 in A1:
            k = nu0 * m1
            i = bisect_right(f2, tn // (td * k))
            small += delta * w1 * k * pwm[i]
            large += delta * w1 * (pw[-1] - pw[i])
    return Fraction(small, T * phi) / l2 + Fraction(large, 2 * M)


def direct_sum(primes: list[int], terms: int = 40) -> Fraction:
    """The bound straight from its definition, by enumerating a mod M (tiny sets only)."""
    M = prod(primes)
    phi = prod(q - 1 for q in primes)
    T = 1
    for q in primes:
        T = lcm(T, ord2(q))
    l2 = log2_lower(terms)
    tot = Fraction(0)
    for a in range(M):
        nu = sum(gcd(a - pow(2, r, M), M) == 1 for r in range(T))
        tot += min(Fraction(1, 2 * M), Fraction(nu, T * phi) / l2)
    return tot


def seed_brute(base: list[int]) -> dict[int, int]:
    M0 = prod(base)
    T0 = 1
    for q in base:
        T0 = lcm(T0, ord2(q))
    h: Counter[int] = Counter()
    for a in range(M0):
        h[sum(gcd(a - pow(2, r, M0), M0) == 1 for r in range(T0))] += 1
    return dict(h)


def self_test() -> None:
    for d, qs in ((11, [23, 89]), (23, [47]), (5, [31]), (3, [7])):
        require(cluster_poly(d, qs) == cluster_poly_brute(d, qs), f"cluster_poly {d}")
    # end to end against the definition: base {3,5} (T0=4), clusters 3:[7], 11:[23]
    base, cl = [3, 5], {3: [7], 11: [23]}
    got = upper_sum(seed_brute(base), base, cl)
    require(got == direct_sum([3, 5, 7, 23]), "upper_sum != direct enumeration")
    for s in range(3):
        require(upper_sum(seed_brute(base), base, cl, split=s) == got, "split changes sum")
    require(Fraction(69, 100) < log2_lower(40) < Fraction(70, 100), "log 2 bound")
    print("verifier self-tests passed", file=sys.stderr)


def check_moments(hist: dict[int, int], base: list[int]) -> None:
    """Independent check of the seed histogram (does not use the DP): its first two
    moments must equal the ones computed straight from the residue conditions,
        sum delta(nu) nu   = T0 prod (q-1),
        sum delta(nu) nu^2 = T0 sum_{k mod T0} prod_q (q - 2 + [ord_q(2) | k]),
    since a pair of exponents r, r+k is blocked by q at 2 residues, or 1 if ord_q | k."""
    T0 = 1
    for q in base:
        T0 = lcm(T0, ord2(q))
    ds = [(q, ord2(q)) for q in base]
    m1 = sum(c * nu for nu, c in hist.items())
    require(m1 == T0 * prod(q - 1 for q in base), "first moment of seed histogram is wrong")
    m2 = sum(c * nu * nu for nu, c in hist.items())
    ref = T0 * sum(prod(q - 2 + (k % d == 0) for q, d in ds) for k in range(T0))
    require(m2 == ref, "second moment of seed histogram is wrong")
    print("seed histogram moments 1 and 2 match the residue conditions", file=sys.stderr)


def load_hist(path: Path) -> dict[int, int]:
    h: dict[int, int] = {}
    for line in path.read_text().splitlines():
        if line.strip():
            nu, c = map(int, line.split())
            require(nu not in h and c > 0, f"bad histogram line {line!r}")
            h[nu] = c
    return h


def dec(f: Fraction, prec: int = 60) -> Decimal:
    getcontext().prec = prec
    return Decimal(f.numerator) / Decimal(f.denominator)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cert", type=Path)
    ap.add_argument("seed", type=Path)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--json-out", type=Path)
    ap.add_argument("--previous", default="0.490180063290061", help="bound to compare with")
    a = ap.parse_args()
    if a.self_test:
        self_test()
    cert = json.loads(a.cert.read_text())
    base = cert["base_order"]
    clusters = {int(d): qs for d, qs in cert["clusters"].items()}
    hist = load_hist(a.seed)
    check_moments(hist, base)
    S = upper_sum(hist, base, clusters)
    allp = sorted(base + [q for qs in clusters.values() for q in qs])
    T = 1
    for q in allp:
        T = lcm(T, ord2(q))
    scale = 10**15
    shown = Fraction(-((-S.numerator * scale) // S.denominator), scale)  # round up
    prev = Fraction(a.previous)
    lines = [
        f"{len(allp)}-prime Romanoff upper-density certificate",
        f"primes = {allp}",
        f"base (processing order) = {base}",
        f"clusters = {json.dumps({str(d): qs for d, qs in clusters.items()})}",
        f"M = {prod(allp)}",
        f"phi(M) = {prod(q - 1 for q in allp)}",
        f"T = {T}",
        f"seed histogram sha256 = {hashlib.sha256(a.seed.read_bytes()).hexdigest()}",
        f"seed histogram entries = {len(hist)}",
        f"rational upper sum = {dec(S, 100)}",
        f"proved: upper density < {dec(shown, 16)}",
        f"improvement over {a.previous} = {dec(prev - S, 30)}",
    ]
    print("\n".join(lines))
    if a.json_out:
        a.json_out.write_text(json.dumps({
            "primes": allp, "base_order": base, "clusters": cert["clusters"], "T": T,
            "upper_sum_numerator": S.numerator, "upper_sum_denominator": S.denominator,
            "upper_sum_decimal": str(dec(S, 100)), "proved_below": str(dec(shown, 16)),
        }, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
