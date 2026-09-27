"""Independent Dn shell enumeration in ambient integer coordinates.

No imports from LDL, the criterion or witness code. Enumerates every integral
v with even coordinate sum and sum v_i^2 <= bound, then inverts the basis.
"""

from collections import Counter, defaultdict
from math import isqrt


def ambient_histograms(rank, bound):
    if type(rank) is not int or rank < 2 or type(bound) is not int or bound < 0:
        raise ValueError("invalid rank or bound")
    histograms = defaultdict(Counter)
    v = []
    leaves = accepted = 0

    def visit(remaining, coordinate_sum):
        nonlocal leaves, accepted
        if len(v) == rank:
            leaves += 1
            if coordinate_sum % 2:
                return
            y, prefix = [], 0
            for j in range(rank - 2):
                prefix += v[j]
                y.append(prefix)
            last = coordinate_sum // 2
            y.extend([last - v[-1], last])
            N = sum(x * x for x in v)
            am = sum((x % 2) << i for i, x in enumerate(y))
            xm = sum(((x // 2) % 2) << i for i, x in enumerate(y))
            histograms[(am, N)][xm] += 1
            accepted += 1
            return
        radius = isqrt(remaining)
        for value in range(-radius, radius + 1):
            v.append(value)
            visit(remaining - value * value, coordinate_sum + value)
            v.pop()

    visit(bound, 0)
    return histograms, {"rank": rank, "norm4_bound": bound,
                        "ambient_leaves": leaves, "lattice_vectors": accepted}
