"""Integer Euclidean enumeration for E6, E7 and E8; no LDL imports.

Simple roots are Bourbaki (1,3,4,...,n,2). Each column of C is twice a
simple root, and C^t C=4G. See E7-E8-MILESTONE.md for completeness.
"""

from collections import Counter, defaultdict
from math import isqrt


def doubled_roots(rank):
    if rank not in (6, 7, 8):
        raise ValueError("rank must be 6, 7 or 8")
    columns = [[1, -1, -1, -1, -1, -1, -1, 1], [-2, 2, 0, 0, 0, 0, 0, 0]]
    for j in range(2, rank - 1):
        column = [0] * 8
        column[j - 1], column[j] = -2, 2
        columns.append(column)
    columns.append([2, 2, 0, 0, 0, 0, 0, 0])
    return [list(row) for row in zip(*columns)]


def gram(rank):
    C = doubled_roots(rank)
    raw = [[sum(row[i] * row[j] for row in C) for j in range(rank)] for i in range(rank)]
    if any(x % 4 for row in raw for x in row):
        raise ValueError("nonintegral root Gram")
    return [[x // 4 for x in row] for row in raw]


def ambient_histograms(rank, bound):
    """Complete counts by (a_mask, norm4, x_mask), in an unshifted sphere."""
    if type(bound) is not int or bound < 0:
        raise ValueError("bound must be nonnegative integral")
    C = doubled_roots(rank)
    weight, budget = 9 - rank, 4 * bound
    histograms = defaultdict(Counter)
    leaves = accepted = 0
    for h in range(-isqrt(budget // weight), isqrt(budget // weight) + 1):
        z = []

        def visit(remaining):
            nonlocal leaves, accepted
            if len(z) == rank - 1:
                leaves += 1
                y = [h] + [0] * (rank - 1)
                for j in range(rank - 2, 1, -1):
                    numerator = z[j] + h + (2 * y[j + 1] if j < rank - 2 else 0)
                    if numerator % 2:
                        return
                    y[j] = numerator // 2
                numerator = z[0] + z[1] + 2 * y[2]
                if numerator % 4:
                    return
                y[-1] = numerator // 4
                numerator = h + 2 * y[-1] - z[0]
                if numerator % 2:
                    return
                y[1] = numerator // 2
                ambient = z + [-h] * (8 - rank) + [h]
                if [sum(c * v for c, v in zip(row, y)) for row in C] != ambient:
                    raise ValueError("inverse-coordinate error")
                scaled = sum(v * v for v in z) + weight * h * h
                if scaled % 4:
                    raise ValueError("nonintegral norm")
                am = sum((v % 2) << j for j, v in enumerate(y))
                xm = sum((((v - v % 2) // 2) % 2) << j for j, v in enumerate(y))
                histograms[(am, scaled // 4)][xm] += 1
                accepted += 1
                return
            radius = isqrt(remaining)
            lo = -radius + ((h + radius) % 2)
            for value in range(lo, radius + 1, 2):
                z.append(value)
                visit(remaining - value * value)
                z.pop()

        visit(budget - weight * h * h)
    return histograms, {"ambient_leaves": leaves, "lattice_vectors": accepted,
                        "rank": rank, "norm4_bound": bound}


def direct_coefficients(histograms, rank):
    """Independent direct integer sign summation for each b; no transform."""
    return {key: [sum(count * (-1 if (xm & bm).bit_count() % 2 else 1)
                      for xm, count in histogram.items()) for bm in range(1 << rank)]
            for key, histogram in histograms.items()}

