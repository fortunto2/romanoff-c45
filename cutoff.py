"""Least histogram cut for a certificate: every nu0 > cut is in the trivial regime for all cluster
atoms, i.e. (cut + 1) * m_min >= thr with thr = T phi(M) lambda / (2M). verify_cert.py checks the same
inequality, so a generator run with --cut equal to this value produces an acceptable aggregate.

usage: python3 cutoff.py CERT.json
"""

import json
import sys
from fractions import Fraction
from math import lcm, prod
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_cert import cluster_poly, log2_lower, ord2


def least_cut(base: list[int], clusters: dict[int, list[int]]) -> int:
    allp = base + [q for qs in clusters.values() for q in qs]
    T = lcm(*(ord2(q) for q in allp))
    thr = Fraction(T * prod(q - 1 for q in allp), 2 * prod(allp)) * log2_lower(40)
    m_min = prod(min(cluster_poly(d, qs)) for d, qs in clusters.items())
    return -(-thr.numerator // (thr.denominator * m_min)) - 1


if __name__ == "__main__":
    c = json.loads(Path(sys.argv[1]).read_text())
    print(least_cut(c["base_order"], {int(d): v for d, v in c["clusters"].items()}))
