"""Uniform Dn criteria, exact witness recipes and compact certificate streams.

See DN-COUNTING.md. Ordinary characteristic coordinates use the simple-root
basis and its dual, while A and C below are ambient vectors 2xi and 4delta.
"""

from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import comb

from exact_theta import (identity, integer, inverse, ldl, matvec, multiply,
                         require, transpose)
from exceptional_batch import bits, ldl_histograms, walsh, histogram_hash


def basis(rank):
    integer(rank, "rank")
    require(rank >= 2, "Dn requires rank >=2")
    B = [[0] * rank for _ in range(rank)]
    for i in range(rank - 1):
        B[i][i], B[i + 1][i] = 1, -1
    B[rank - 2][-1] = B[rank - 1][-1] = 1
    return B


def gram(rank):
    B = basis(rank)
    return multiply(transpose(B), B)


def formula(rank):
    integer(rank, "rank")
    require(rank >= 2, "rank must be >=2")
    return (4 ** (rank - 1) - 3 ** rank + 3 * 2 ** (rank - 1) - 1
            + (comb(rank, rank // 2) // 2 if rank % 2 == 0 else 0))


def product_zero_count(rank):
    return sum(comb(rank, s) * (2 ** (s - 1) - 2) * 2 ** (rank - s)
               for s in range(2, rank + 1, 2))


def ambient_a(am, rank):
    return matvec(basis(rank), bits(am, rank))


def ambient_c(bm, rank):
    """C=2B_ambient=4delta; B_basis^t C=2 b_dual."""
    b = bits(bm, rank)
    C = [0] * rank
    C[-1] = b[-1] - b[-2]
    running = 0
    for i in range(rank - 2, -1, -1):
        running += b[i]
        C[i] = C[-1] + 2 * running
    return C


def characteristic_tables(rank):
    B = basis(rank)
    primals, duals = [], []
    for am in range(1 << rank):
        A = matvec(B, bits(am, rank))
        S = sum((x % 2) << i for i, x in enumerate(A))
        e = ((sum(A) - S.bit_count()) // 2) % 2
        primals.append((A, S, e))
    for bm in range(1 << rank):
        C = ambient_c(bm, rank)
        h = C[0] % 2
        require(all(x % 2 == h for x in C), "dual vector has mixed fractional parts")
        require(matvec(transpose(B), C) == [2 * x for x in bits(bm, rank)], "dual-coordinate error")
        W = sum((((x - h) // 2) % 2) << i for i, x in enumerate(C))
        duals.append((C, h, W))
    return primals, duals


def zero_recipe(am, bm, rank, tables):
    """Return a witness recipe for an even zero, or None for a nonzero class."""
    require((am & bm).bit_count() % 2 == 0, "even pair required")
    _, S, e = tables[0][am]
    _, h, W = tables[1][bm]
    if h:
        return None
    s, r = S.bit_count(), (S & W).bit_count()
    if s and 0 < r < s:
        i = (S & W & -(S & W)).bit_length() - 1
        other = S & ~W
        j = (other & -other).bit_length() - 1
        return [1, i, j]  # flip these two ambient coordinates
    if s == 0 and e == 1 and 2 * W.bit_count() == rank:
        return [2]  # interchange W with its complement
    return None


def signed_permutation(recipe, rank, W):
    require(isinstance(recipe, list) and recipe, "witness recipe required")
    kind = integer(recipe[0], "recipe kind")
    permutation, signs = list(range(rank)), [1] * rank
    if kind == 1:
        require(len(recipe) == 3, "flip recipe needs two coordinates")
        i, j = integer(recipe[1], "i"), integer(recipe[2], "j")
        require(0 <= i < rank and 0 <= j < rank and i != j, "invalid flip coordinates")
        signs[i] = signs[j] = -1
    elif kind == 2:
        require(len(recipe) == 1, "swap recipe has no extra parameters")
        left = [i for i in range(rank) if (W >> i) & 1]
        right = [i for i in range(rank) if not ((W >> i) & 1)]
        require(len(left) == len(right), "balanced swap requires equal sizes")
        for i, j in zip(left, right):
            permutation[i], permutation[j] = j, i
    else:
        require(False, "unknown witness recipe")
    require(sorted(permutation) == list(range(rank)), "not a permutation")
    return permutation, signs


def verify_ambient_witness(A, C, permutation, signs):
    """Check stabilizers and epsilon in integers; no matrix inversion per class."""
    SA = [s * A[j] for s, j in zip(signs, permutation)]
    SC = [s * C[j] for s, j in zip(signs, permutation)]
    da = [x - y for x, y in zip(SA, A)]
    dc = [x - y for x, y in zip(SC, C)]
    require(all(x % 2 == 0 for x in da) and sum(da) % 4 == 0,
            "witness does not stabilize xi modulo Dn")
    residues = {x % 4 for x in dc}
    require(residues in ({0}, {2}), "witness does not stabilize delta modulo Dn*")
    numerator = sum(x * y for x, y in zip(A, dc))
    require(numerator % 4 == 0 and (numerator // 4) % 2 == 1, "epsilon is not one")


def header(rank):
    bound = max(4, 2 * (rank // 2))
    return {"schema": "work7-dn-stream-v1", "rank": rank, "G": gram(rank),
            "representatives": "binary-primal-and-dual-masks", "norm4_bound": bound,
            "record_format": "[a,b,0,norm4,coefficient] or [a,b,1,i,j] or [a,b,2]"}


def validate_header(document):
    require(isinstance(document, dict) and document.get("schema") == "work7-dn-stream-v1", "unknown schema")
    rank = integer(document.get("rank"), "rank")
    require(2 <= rank <= 10, "finite audit domain is D2 through D10")
    require(document.get("G") == gram(rank), "wrong Dn Gram matrix")
    # Do not accept bools or floats that happen to compare equal to integers.
    require(all(type(x) is int for row in document["G"] for x in row), "Gram entries must be integers")
    require(document.get("representatives") == "binary-primal-and-dual-masks", "wrong representatives")
    bound = integer(document.get("norm4_bound"), "norm4_bound")
    require(bound == max(4, 2 * (rank // 2)), "unexpected finite search bound")
    _, D = ldl(document["G"])
    det = Fraction(1)
    for d in D:
        det *= d
    require(det == 4, "Dn determinant must be four")
    return rank, bound


def expected_pairs(rank):
    for am in range(1 << rank):
        for bm in range(1 << rank):
            if (am & bm).bit_count() % 2 == 0:
                yield am, bm


def coefficient_tables(document, backend="ldl", max_nodes=1_000_000):
    rank, bound = validate_header(document)
    if backend == "ldl":
        histograms, stats = ldl_histograms(document["G"], bound, max_nodes)
    elif backend == "ambient":
        from dn_ambient import ambient_histograms
        histograms, stats = ambient_histograms(rank, bound)
    else:
        require(False, "unknown backend")
    coefficients = {key: walsh(histogram, rank) for key, histogram in histograms.items()}
    # All possible b are evaluated together, including the odd parity check.
    for (am, _), values in coefficients.items():
        require(all(c == 0 for bm, c in enumerate(values) if (am & bm).bit_count() % 2),
                "odd characteristic has a nonzero coefficient")
    norms = [[] for _ in range(1 << rank)]
    for am, N in sorted(coefficients):
        norms[am].append(N)
    return coefficients, norms, histograms, stats


def records(document, coefficients, norms):
    rank = document["rank"]
    tables = characteristic_tables(rank)
    for am, bm in expected_pairs(rank):
        recipe = zero_recipe(am, bm, rank, tables)
        if recipe is not None:
            yield [am, bm, *recipe]
        else:
            N = next((N for N in norms[am] if coefficients[(am, N)][bm]), None)
            require(N is not None, f"nonzero class unresolved: a={am}, b={bm}")
            yield [am, bm, 0, N, coefficients[(am, N)][bm]]


def verify_records(document, rows, coefficients, norms, histograms, stats, backend):
    rank = document["rank"]
    tables = characteristic_tables(rank)
    iterator = iter(rows)
    digest = hashlib.sha256()
    counts, witness_counts = Counter(), Counter()
    permutation_cache = {}
    count = 0
    for am, bm in expected_pairs(rank):
        row = next(iterator, None)
        require(isinstance(row, list) and all(type(x) is int for x in row), "integral record required")
        require(len(row) >= 3 and row[:2] == [am, bm], "missing, duplicate or out-of-order class")
        kind = row[2]
        if kind == 0:
            require(len(row) == 5, "shell record must have five entries")
            N, c = row[3:]
            require(0 <= N <= document["norm4_bound"] and c != 0, "nonzero shell within range required")
            require(coefficients.get((am, N), [0] * (1 << rank))[bm] == c, "coefficient mismatch")
            require(zero_recipe(am, bm, rank, tables) is None, "zero criterion contradicts nonzero evidence")
            counts[0] += 1
            witness_counts[N] += 1
        else:
            require(kind in (1, 2), "unknown record kind")
            require(row[2:] == zero_recipe(am, bm, rank, tables), "recipe does not match the uniform criterion")
            A = tables[0][am][0]
            C, _, W = tables[1][bm]
            cache_key = tuple(row[2:]) + ((W,) if kind == 2 else ())
            if cache_key not in permutation_cache:
                permutation_cache[cache_key] = signed_permutation(row[2:], rank, W)
            permutation, signs = permutation_cache[cache_key]
            verify_ambient_witness(A, C, permutation, signs)
            require(all(coefficients[(am, N)][bm] == 0 for N in norms[am]), "witness contradicts a shell")
            counts[kind] += 1
        digest.update((json.dumps(row, separators=(",", ":")) + "\n").encode("ascii"))
        count += 1
    require(next(iterator, None) is None, "extra record after complete coverage")
    zeros = counts[1] + counts[2]
    require(zeros == formula(rank), "zero count disagrees with closed formula")
    require(counts[1] == product_zero_count(rank), "product-zero combinatorial count mismatch")
    require(counts[2] == (comb(rank, rank // 2) // 2 if rank % 2 == 0 else 0), "balanced count mismatch")
    return {"schema": "work7-dn-result-v1", "verdict": "classification_verified", "rank": rank,
            "backend": backend, "determinant": 4, "coverage_complete": True,
            "total_characteristics": 4 ** rank, "even_characteristics": count,
            "odd_characteristics_vanish_by_parity": 4 ** rank - count,
            "even_nonvanishing": counts[0], "even_vanishing": zeros,
            "product_zero_classes": counts[1], "balanced_cancellation_classes": counts[2],
            "formula_value": formula(rank), "unresolved": 0,
            "maximum_nonzero_witness_norm": str(Fraction(max(witness_counts), 4)),
            "witness_norm4_counts": {str(N): witness_counts[N] for N in sorted(witness_counts)},
            "records_sha256": digest.hexdigest(), "shell_histogram_sha256": histogram_hash(histograms),
            "enumeration": stats, "full_automorphism_group_claimed": False}


def expand_record(document, row):
    rank = document["rank"]
    am, bm, kind = row[:3]
    if kind == 0:
        evidence = {"kind": "shell", "norm4": row[3], "coefficient": row[4]}
    else:
        C = ambient_c(bm, rank)
        h = C[0] % 2
        W = sum((((x - h) // 2) % 2) << i for i, x in enumerate(C))
        permutation, signs = signed_permutation(row[2:], rank, W)
        P = [[signs[i] if j == permutation[i] else 0 for j in range(rank)] for i in range(rank)]
        B = basis(rank)
        raw = multiply(multiply(inverse(B), P), B)
        require(all(x.denominator == 1 for line in raw for x in line), "nonintegral lattice action")
        evidence = {"kind": "symmetry", "T": [[int(x) for x in line] for line in raw], "order": 2}
    return {"schema": "work7-theta-v1", "name": f"D{rank}: a={am}, b={bm}", "G": document["G"],
            "a": bits(am, rank), "b": bits(bm, rank), "evidence": evidence}
