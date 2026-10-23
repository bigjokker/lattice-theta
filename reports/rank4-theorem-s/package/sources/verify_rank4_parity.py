"""Independent replay of the rank-four parity results (RANK4-PARITY-CONJECTURE.md).

Nothing here imports rank4_parity_search.  Candidate domains are regenerated
with a loose-box diagonal scanner, groups and transports with the box
isometry search of verify_rank4_feasibility, coefficients with the box
enumerations of verify_lattice_census (small bounds) and verify_rank4_level4
(large cutoffs), zeros with exact_theta.symmetry_data or the modular cutoff
replay.  A node limit raises; nothing partial is certified.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import floor, prod
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, dot, identity, integer,
                         inverse, ldl, matrix, matvec, multiply, require, symmetry_data,
                         transpose, vector)
from verify_lattice_census import box_shells, box_vectors
from verify_rank4_feasibility import box_isometries, cutoff_replay, det4, modular_check
from verify_rank4_level4 import box_signed_shells

BITS = [list(v) for v in product((0, 1), repeat=4)]
EVEN = [(a, b) for a in BITS if any(a) for b in BITS if any(b) and dot(a, b) % 2 == 0]
I4 = identity(4)
MINUS_I4 = [[-x for x in row] for row in I4]
DET_FILTERS = {"all": lambda d: True, "odd": lambda d: d % 2 == 1,
               "divisible_by_4": lambda d: d % 4 == 0, "not_divisible_by_4": lambda d: d % 4 != 0}
BIG_BOUND = 2000


def even(G):
    return all(G[i][i] % 2 == 0 for i in range(4))


FAST_BOX = [False]


def shells(G, a, b, bound, max_nodes):
    """Box backends: the rectangular box of verify_lattice_census for small bounds, the outer box with
    exact inner solve of verify_rank4_level4 for large bounds (or everywhere with --fast-box)."""
    if bound <= BIG_BOUND and not FAST_BOX[0]:
        return box_shells(G, a, b, bound, max_nodes)
    return box_signed_shells(G, a, b, bound, max_nodes)[0]


# ----------------------------------------------------------------------------
# independent candidate domain
# ----------------------------------------------------------------------------

def independent_candidates(D0, even_only, det_filter, dmin=1):
    """Loose cross-entry boxes and direct diagonal scanning (RANK4-FEASIBILITY.md section 1)."""
    integer(D0, "D0")
    require(1 <= dmin <= D0, "wrong domain")
    det_ok = DET_FILTERS[det_filter]
    ratio = Fraction(3, 4)
    upper = [Fraction(D0) / ratio ** (6 - i) for i in range(4)]
    found = []

    def extend(G):
        k = len(G)
        mu, d = ldl(G)
        inv = inverse(G)
        prefix = prod(d)
        radii = [floor(sum(abs(mu[j][l]) * d[l] for l in range(j + 1)) / 2) for j in range(k)]
        for cross in product(*(range(-r, r + 1) for r in radii)):
            coeff = matvec(inv, cross)
            projected = [sum(coeff[i] * mu[i][j] for i in range(j, k)) for j in range(k)]
            if any(abs(x) > Fraction(1, 2) for x in projected):
                continue
            q = dot(cross, coeff)
            for diagonal in range(G[0][0], floor(upper[k] + sum(upper[:k]) / 4) + 1):
                if even_only and diagonal % 2:
                    continue
                pivot = diagonal - q
                if pivot <= 0 or pivot < ratio * d[-1]:
                    continue
                r = 4 - k
                if prefix * pivot ** r * ratio ** (r * (r - 1) // 2) > D0:
                    break
                row = [list(old) + [x] for old, x in zip(G, cross)] + [list(cross) + [diagonal]]
                if k == 3:
                    delta = prefix * pivot
                    if delta.denominator == 1 and dmin <= delta <= D0 and det_ok(int(delta)):
                        found.append(row)
                else:
                    extend(row)

    for A in range(1, D0 + 1):
        if even_only and A % 2:
            continue
        if ratio ** 6 * A ** 4 <= D0:
            extend([[A]])
    return sorted(found)


# ----------------------------------------------------------------------------
# small independent helpers
# ----------------------------------------------------------------------------

def splitting_involution(group):
    for T in group:
        if T != I4 and T != MINUS_I4 and multiply(T, T) == I4 and all((x + y) % 2 == 0 for rT, rI in zip(T, I4) for x, y in zip(rT, rI)):
            return T
    return None


def indecomposable_components(G, max_nodes):
    """Eichler components over the box sphere of radius max G_ii (see the builder docstring).

    Same vectors, same two criteria as indecomposable_components_reference (kept for tests): a vector v is
    indecomposable unless some u with N(u) < N(v) has B(u, v) = N(u); components join indecomposable
    vectors with nonzero inner product.  G v is computed once per vector, the reducibility scan stops at
    the first u with N(u) >= N(v) (vectors sorted by norm), and components come from a union-find.
    """
    R = max(G[i][i] for i in range(4))
    vecs = [tuple(v) for v in box_vectors(G, R, max_nodes) if any(v)]
    gv = {v: tuple(matvec(G, v)) for v in vecs}
    norm = {v: dot(v, gv[v]) for v in vecs}
    ordered = sorted(vecs, key=norm.get)
    indec = []
    for v in vecs:
        if v < tuple(-x for x in v):
            continue
        nv, g = norm[v], gv[v]
        reducible = False
        for u in ordered:
            nu = norm[u]
            if nu >= nv:
                break
            if u[0] * g[0] + u[1] * g[1] + u[2] * g[2] + u[3] * g[3] == nu:
                reducible = True
                break
        if not reducible:
            indec.append(v)
    parent = list(range(len(indec)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, v in enumerate(indec):
        g = gv[v]
        for j in range(i):
            u = indec[j]
            if u[0] * g[0] + u[1] * g[1] + u[2] * g[2] + u[3] * g[3]:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    return len({find(i) for i in range(len(indec))})


def indecomposable_components_reference(G, max_nodes):
    """The original implementation, kept as the regression oracle for indecomposable_components."""
    R = max(G[i][i] for i in range(4))
    vecs = [tuple(v) for v in box_vectors(G, R, max_nodes) if any(v)]
    norm = {v: dot(v, matvec(G, v)) for v in vecs}
    indec = []
    for v in vecs:
        if v < tuple(-x for x in v):
            continue
        Gv = matvec(G, v)
        if not any(norm[u] < norm[v] and dot(u, Gv) == norm[u] for u in vecs):
            indec.append(v)
    labels = list(range(len(indec)))
    changed = True
    while changed:
        changed = False
        for i in range(len(indec)):
            for j in range(len(indec)):
                if labels[i] != labels[j] and dot(indec[i], matvec(G, indec[j])) != 0:
                    m = min(labels[i], labels[j])
                    labels[i] = labels[j] = m
                    changed = True
    return len(set(labels))


def q2(G, x):
    n = dot(x, matvec(G, x))
    require(n % 2 == 0, "odd norm in an even lattice")
    return (n // 2) % 2


def f2_check(G, a, b, stored):
    """Recompute the F2 geometry by brute force over the 16 classes."""
    rad = [x for x in BITS if all(dot(x, matvec(G, y)) % 2 == 0 for y in BITS)]
    cs = [x for x in BITS if all((u - v) % 2 == 0 for u, v in zip(matvec(G, x), b))]
    singular = sum(1 for x in BITS if any(x) and q2(G, x) == 0)
    rank = 4 - (len(rad).bit_length() - 1)
    expect = {"rank_mod_2": rank, "radical_dimension": 4 - rank, "nonzero_singular_vectors": singular,
              "a_in_radical": list(a) in rad, "Q_a": q2(G, a), "b_in_image_mod_2": bool(cs),
              "radical_totally_singular": all(q2(G, x) == 0 for x in rad), "pairing_a_b_mod_2": dot(a, b) % 2}
    if rank == 4:
        expect["type"] = "plus" if singular == 9 else "minus"
        require(singular in (5, 9), "nondegenerate quadratic space with a wrong singular count")
    else:
        expect["type"] = "degenerate"
    if cs:
        expect["Q_c_values"] = sorted({q2(G, c) for c in cs})
        expect["a_equals_some_c"] = list(a) in cs
        expect["span_a_c_totally_singular"] = expect["Q_a"] == 0 and expect["Q_c_values"] == [0] and list(a) not in cs
    for key, value in expect.items():
        require(stored.get(key) == value, "F2 geometry mismatch: " + key)
    require(sorted(stored["c_solutions"]) == sorted(cs) and sorted(stored["radical"]) == sorted(rad), "F2 sets mismatch")


def stabilizer(G, a, b, group, duals):
    indices, signs = [], []
    for i, (T, dual) in enumerate(zip(group, duals)):
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        transported = matvec(dual, b)
        if any((x - y) % 2 for x, y in zip(transported, b)):
            continue
        shift = [(x - y) / 2 for x, y in zip(transported, b)]
        require(all(x.denominator == 1 for x in shift), "nonintegral shift")
        indices.append(i)
        signs.append(int(dot(a, shift)) % 2)
    return indices, signs


def check_kernel(G, a, b, kernel, max_nodes):
    """Sign kernel, classes, coset lattices and level of a cutoff zero."""
    K = matrix(kernel["kernel_basis"], "kernel basis", 4)
    require(abs(det4(K)) == 2 and all(dot(col, b) % 2 == 0 for col in transpose(K)), "wrong kernel basis")
    G0 = multiply(multiply(transpose(K), G), K)
    require(kernel["G0"] == G0 and kernel["d0"] == det4(G0) == 4 * det4(G), "wrong kernel Gram")
    t = kernel["complement_t"]
    require(dot(t, b) % 2 == 1, "complement not odd")
    Kinv = inverse(K)
    w1, w2 = matvec(Kinv, a), matvec(Kinv, [p + 2 * q for p, q in zip(a, t)])
    require(all(x.denominator == 1 for x in w1 + w2), "classes not in the kernel")
    v1, v2 = [int(x) % 2 for x in w1], [int(x) % 2 for x in w2]
    require(kernel["classes"] == {"v1": v1, "v2": v2}, "wrong classes")
    group0 = box_isometries(G0, G0, max_nodes)
    require(kernel["aut_L0_order"] == len(group0), "wrong Aut(L0) order")

    def cls(v):
        return tuple(x % 2 for x in v)
    mapping = [T for T in group0 if cls(matvec(T, v1)) == tuple(v2)]
    swapping = [T for T in mapping if cls(matvec(T, v2)) == tuple(v1)]
    require(kernel["elements_mapping_v1_to_v2"] == len(mapping) and kernel["elements_swapping_v1_v2"] == len(swapping),
            "wrong Aut(L0) action")
    for name, v in (("O_v1", v1), ("O_v2", v2)):
        c = kernel["coset_lattices"][name]
        C = matrix(c["doubled_basis_in_L0"], "coset basis", 4)
        require(abs(det4(C)) == 8 and all(all((x - y) % 2 == 0 for x, y in zip(col, [0] * 4)) or cls(col) == tuple(v) for col in transpose(C)),
                "wrong coset basis")
        require(c["gram_times_4"] == multiply(multiply(transpose(C), G0), C), "wrong coset Gram")
        red = c["reduced"]
        Bm = matrix(red["basis"], "reduced basis", 4)
        U = multiply(inverse(C), Bm)
        require(all(x.denominator == 1 for row in U for x in row) and abs(det4([[int(x) for x in row] for row in U])) == 1, "reduced basis not a basis")
        require(red["gram"] == multiply(multiply(transpose(Bm), G0), Bm) and red["determinant"] == det4(red["gram"]), "wrong reduced Gram")
    if swapping:
        level = 1
    elif mapping:
        level = 2
    else:
        iso = box_isometries(kernel["coset_lattices"]["O_v1"]["reduced"]["gram"],
                             kernel["coset_lattices"]["O_v2"]["reduced"]["gram"], max_nodes)
        level = 3 if iso else 4
    require(kernel["level"] == level, "wrong level")
    return level


def check_decisions(G, group, par, max_nodes, enumeration_nodes, external_ok=True):
    """Replay a decide_parity document for the lattice G with the replayed group."""
    require(par["G"] == G and par["aut_order"] == len(group), "parity document for another lattice")
    duals = [transpose(inverse(T)) for T in group]
    decisions = par["decisions"]
    require([(d["a"], d["b"]) for d in decisions] == EVEN, "characteristic coverage")
    verdicts = Counter()
    zeros, firsts, minimal_decides = 0, [], 0
    levels = Counter()
    for d in decisions:
        a, b = vector(d["a"], "a", 4), vector(d["b"], "b", 4)
        indices, signs = stabilizer(G, a, b, group, duals)
        require(d["stabilizer"] == {"indices": indices, "epsilon": signs}, "wrong stabilizer")
        if even(G) and "f2" in d:
            f2_check(G, a, b, d["f2"])
        verdict = d["verdict"]
        verdicts[verdict] += 1
        if verdict == "zero_by_symmetry":
            T = matrix(d["witness"], "witness", 4)
            require(T in group and 1 in signs and symmetry_data(G, a, b, T)["epsilon"] == 1, "bad symmetry witness")
            require(d["level"] == 1, "symmetry zero must be level 1")
            zeros += 1
            levels[1] += 1
        elif verdict == "nonzero":
            require(1 not in signs, "nonzero despite a symmetry witness")
            bound = integer(d["norm4_bound"], "bound")
            got = shells(G, a, b, bound, enumeration_nodes)
            require(got == d["shells"], "coefficient histogram mismatch")
            first = next(s for s in got if s["signed_coefficient"])
            require(d["first_nonzero"] == {"norm4": first["norm4"], "coefficient": first["signed_coefficient"]}, "wrong first nonzero")
            m = next(s for s in got if s["vector_count"])
            require(d["minimal_shell"] == {"norm4": m["norm4"], "vector_count": m["vector_count"], "signed_coefficient": m["signed_coefficient"]}, "wrong minimal shell")
            firsts.append(first["norm4"])
            minimal_decides += first["norm4"] == m["norm4"]
            if "cutoff_certificate" in d:
                cert = d["cutoff_certificate"]
                require(cert["G"] == G and cert["a"] == a and cert["b"] == b, "cutoff for another characteristic")
                checked = cutoff_replay(cert, enumeration_nodes)
                require(checked["certificate_replay_complete"] and checked["verdict"] == "proved_nonzero", "cutoff nonzero replay")
        elif verdict == "zero_by_cutoff":
            require(1 not in signs, "cutoff zero despite a symmetry witness")
            cert = d["cutoff_certificate"]
            require(cert["G"] == G and cert["a"] == a and cert["b"] == b, "cutoff for another characteristic")
            B = modular_check(G, a, b, cert["modular"])
            require(cert["run"]["norm4_bound"] >= B and cert["evidence"]["enumeration_complete"] is True, "cutoff not reached")
            got = shells(G, a, b, cert["run"]["norm4_bound"], enumeration_nodes)
            require(got == cert["evidence"]["enumeration"]["shells"], "cutoff histogram mismatch")
            require(all(s["signed_coefficient"] == 0 for s in got), "a nonzero coefficient inside a claimed zero")
            require(cert["evidence"]["verdict"] == "proved_zero_by_modular_cutoff", "wrong cutoff verdict")
            zeros += 1
            if "kernel" in d:
                require(d["level"] == check_kernel(G, a, b, d["kernel"], max_nodes), "level")
                levels[d["level"]] += 1
        elif verdict == "zero_certified_externally":
            require(external_ok and 1 not in signs, "unexpected external zero")
            ext = d["external_certificate"]
            raw = Path(ext["path"]).read_bytes()
            require(hashlib.sha256(raw).hexdigest() == ext["sha256"], "external certificate hash")
            doc = json.loads(raw)
            require(any(c["G"] == G and c["a"] == a and c["b"] == b and c["level"] == ext["level"] for c in doc["certificates"]),
                    "external certificate does not contain this zero")
            zeros += 1
            levels[ext["level"]] += 1
        else:
            require(verdict == "unresolved" and 1 not in signs, "unknown verdict")
    status = "holds" if verdicts["nonzero"] == 105 else "fails" if zeros else "unresolved"
    require(par["P"] == status and par["verdicts"] == dict(verdicts) and par["even_zeros"] == zeros, "wrong P summary")
    require({str(k): v for k, v in par["first_nonzero_norm4_distribution"].items()}
            == {str(k): v for k, v in sorted(Counter(firsts).items())}, "wrong first-nonzero table")
    require(par["largest_first_nonzero_norm4"] == (max(firsts) if firsts else None), "wrong largest first nonzero")
    require(par["minimal_shell_decides"] == minimal_decides, "wrong minimal-shell count")
    if par.get("zero_levels") is not None and "zero_certified_externally" not in verdicts:
        require({int(k): v for k, v in par["zero_levels"].items()} == dict(levels), "wrong level tally")
    return status


# ----------------------------------------------------------------------------
# shared candidate domain for an interval of determinants
# (RANK4-NIGHT-REVIEW-20261022.md section 1; Grok Bot's per-determinant re-walk)
# ----------------------------------------------------------------------------

RATIO = Fraction(3, 4)


@lru_cache(maxsize=200_000)
def _edge(prefix, cross):
    """The D0-independent part of one step of independent_candidates.extend, for a fixed prefix and cross row."""
    G = [list(r) for r in prefix]
    k = len(G)
    mu, d = ldl(G)
    coeff = matvec(inverse(G), list(cross))
    radii = [floor(sum(abs(mu[j][l]) * d[l] for l in range(j + 1)) / 2) for j in range(k)]
    projected = [sum(coeff[i] * mu[i][j] for i in range(j, k)) for j in range(k)]
    valid = all(abs(x) <= r for x, r in zip(cross, radii)) and all(abs(x) <= Fraction(1, 2) for x in projected)
    return valid, dot(list(cross), coeff), prod(d), d[-1]


def admissible(row, D0, even_only, det_filter, dmin):
    """True iff independent_candidates(D0, even_only, det_filter, dmin) emits this row.

    The scanner builds a form row by row, so a row has exactly one construction path (its leading
    principal blocks).  Every test of extend() is re-applied along that path with this D0: the start
    condition on A, the cross-entry box and the |projected| <= 1/2 test, the diagonal range and parity,
    the pivot conditions, the D0 break (monotone in the pivot, hence a test at this diagonal), and the
    final determinant window and filter.
    """
    det_ok = DET_FILTERS[det_filter]
    upper = [Fraction(D0) / RATIO ** (6 - i) for i in range(4)]
    A = row[0][0]
    if not (1 <= A <= D0) or (even_only and A % 2) or RATIO ** 6 * A ** 4 > D0:
        return False
    for k in range(1, 4):
        ok, q, prefix, last = _edge(tuple(tuple(r[:k]) for r in row[:k]), tuple(row[k][:k]))
        diagonal = row[k][k]
        if not ok or (even_only and diagonal % 2) or not A <= diagonal <= floor(upper[k] + sum(upper[:k]) / 4):
            return False
        pivot = diagonal - q
        r = 4 - k
        if pivot <= 0 or pivot < RATIO * last or prefix * pivot ** r * RATIO ** (r * (r - 1) // 2) > D0:
            return False
        if k == 3:
            delta = prefix * pivot
            return delta.denominator == 1 and dmin <= delta <= D0 and det_ok(int(delta))
    return False


class SharedDomain:
    """Per-determinant candidate domains from ONE run of independent_candidates over an interval.

    build() scans [min(dets), max(dets)] once, buckets every emitted form by its exact determinant, and
    keeps a form in bucket d only if admissible(form, d, ...) holds, i.e. only if the one-determinant
    scan independent_candidates(d, even_only, det_filter, d) would emit it.  Buckets are produced only
    by this class from this module's own scanner; replay_class_list still compares the bucket with the
    stored candidate list, transports and retention order exactly as before.  Use it for a batch: a
    single determinant is cheaper with the ordinary scan.

    Provenance: this is a TRUSTED INTERNAL object.  The type test in replay_class_list is an API guard
    against passing a plain list, not a proof of provenance: the constructor itself accepts arbitrary
    buckets.  Only domains built by SharedDomain.build in the verifying process should be used; it does not
    validate externally supplied or cached buckets.
    """

    def __init__(self, dets, even_only, det_filter, buckets, info):
        self.dets, self.even_only, self.det_filter = frozenset(dets), even_only, det_filter
        self._buckets, self.info = buckets, info

    @classmethod
    def build(cls, dets, even_only, det_filter):
        dets = sorted(set(dets))
        require(bool(dets) and all(isinstance(d, int) and d >= 1 for d in dets), "shared domain: positive determinants")
        require(det_filter in DET_FILTERS, "shared domain: unknown determinant filter")
        lo, hi = dets[0], dets[-1]
        start = time.perf_counter()
        forms = independent_candidates(hi, even_only, det_filter, lo)
        scan = time.perf_counter() - start
        want = set(dets)
        buckets = {d: [] for d in dets}
        discarded = 0
        for G in forms:
            d = det4(G)
            if d not in want:
                continue
            if admissible(G, d, even_only, det_filter, d):
                buckets[d].append(G)
            else:
                discarded += 1
        info = {"interval": [lo, hi], "determinants": len(dets), "even_only": even_only, "det_filter": det_filter,
                "forms_in_interval_scan": len(forms), "forms_in_buckets": sum(len(b) for b in buckets.values()),
                "discarded_by_per_determinant_guard": discarded, "scan_seconds": round(scan, 2),
                "seconds": round(time.perf_counter() - start, 2),
                "scanner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        return cls(dets, even_only, det_filter, buckets, info)

    def restrict(self, d):
        """The domain for one determinant (for handing to a worker process)."""
        require(d in self.dets, "shared domain: determinant not in this domain")
        return SharedDomain([d], self.even_only, self.det_filter, {d: self._buckets[d]}, dict(self.info, restricted_to=d))

    def candidates(self, D0, even_only, det_filter, dmin):
        require(D0 == dmin and D0 in self.dets, "shared domain: request is not one determinant of this domain")
        require(even_only == self.even_only and det_filter == self.det_filter, "shared domain: flags differ")
        return [list(map(list, G)) for G in self._buckets[D0]]


# ----------------------------------------------------------------------------
# distinctness prefilter: shell counts prove NON-isometry only
# ----------------------------------------------------------------------------

PROFILE_NORM = 16


def shell_profile(G, max_nodes):
    """Number of vectors of each norm 0..PROFILE_NORM, from the complete box of verify_lattice_census.box_vectors.
    An isometry preserves these counts, so different profiles prove two lattices non-isometric."""
    counts = Counter(dot(v, matvec(G, v)) for v in box_vectors(G, PROFILE_NORM, max_nodes))
    return tuple(counts.get(q, 0) for q in range(PROFILE_NORM + 1))


# ----------------------------------------------------------------------------
# class lists and decided class lists (Tasks B, C, E)
# ----------------------------------------------------------------------------

def replay_class_list(cdoc, max_nodes, domain=None):
    """Independent replay of a class list.  domain: an optional SharedDomain that supplies this
    determinant's candidate domain from one interval scan (identical to the per-determinant scan)."""
    require(cdoc["schema"] == "work7-rank4-parity-v1-class-list", "wrong class-list schema")
    D0, dmin = integer(cdoc["determinant_bound"], "D0"), integer(cdoc["determinant_minimum"], "dmin")
    if domain is None:
        forms = independent_candidates(D0, cdoc["even_only"], cdoc["determinant_filter"], dmin)
    else:
        require(isinstance(domain, SharedDomain), "replay: domain must be a SharedDomain")
        forms = domain.candidates(D0, cdoc["even_only"], cdoc["determinant_filter"], dmin)
    if cdoc["odd_lattices_only"]:
        forms = [G for G in forms if not even(G)]
    require(cdoc["candidate_bases"] == forms and cdoc["candidates"] == len(forms), "candidate domain mismatch")
    reps = cdoc["representatives"]
    require(cdoc["classes"] == len(reps), "class count")
    for r in reps:
        G = matrix(r["G"], "representative", 4)
        ldl(G)
        require(r["determinant"] == det4(G) and dmin <= det4(G) <= D0 and DET_FILTERS[cdoc["determinant_filter"]](det4(G)), "representative determinant")
        require(not cdoc["even_only"] or even(G), "parity of representative")
        require(not cdoc["odd_lattices_only"] or not even(G), "parity of representative")
    seen = []
    require(len(cdoc["generation"]) == len(forms), "generation coverage")
    for n, item in enumerate(cdoc["generation"]):
        require(item["candidate"] == n, "generation order")
        idx = integer(item["representative"], "representative index")
        require(0 <= idx < len(reps), "representative index")
        T = matrix(item["T"], "transport", 4)
        require(abs(det4(T)) == 1 and multiply(multiply(transpose(T), reps[idx]["G"]), T) == forms[n], "invalid transport")
        if idx not in seen:
            require(idx == len(seen) and reps[idx]["G"] == forms[n], "retention order")
            seen.append(idx)
    require(seen == list(range(len(reps))), "unused representative")
    searches = separated = full = 0
    profiles = {}
    for i, r in enumerate(reps):
        for s in reps[:i]:
            if r["determinant"] == s["determinant"]:
                searches += 1                 # pairs proved distinct (same meaning as before)
                for x in (r, s):
                    key = json.dumps(x["G"])
                    if key not in profiles:
                        profiles[key] = shell_profile(x["G"], max_nodes)
                if profiles[json.dumps(r["G"])] != profiles[json.dumps(s["G"])]:
                    separated += 1            # different shell counts: proved non-isometric
                    continue
                full += 1                     # equal counts are never taken as isometry: complete search
                require(not box_isometries(r["G"], s["G"], max_nodes), "duplicate isometry class")
    return {"candidates": len(forms), "classes": len(reps), "pairwise_searches": searches,
            "separated_by_shell_counts": separated, "full_isometry_searches": full}


def replay_decided(doc, max_nodes, enumeration_nodes, progress=None):
    require(doc["schema"] == "work7-rank4-parity-v1-decided-classes", "wrong schema")
    start = time.perf_counter()
    listing = replay_class_list(doc["classes_doc"], max_nodes)
    reps = doc["classes_doc"]["representatives"]
    rows = doc["lattices"]
    require([r["G"] for r in rows] == [r["G"] for r in reps], "decided rows do not match the class list")
    holds, fails, unresolved, indec = 0, 0, 0, 0
    firsts, levels = Counter(), Counter()
    failing = []
    for i, row in enumerate(rows):
        G = row["G"]
        group = box_isometries(G, G, max_nodes)
        require(row["automorphisms"] == group and row["aut_order"] == len(group), "full group mismatch")
        require(row["determinant"] == det4(G) and row["even"] == even(G), "row metadata")
        require(row["has_square_root_of_minus_identity"] == any(multiply(T, T) == MINUS_I4 for T in group), "sqrt(-I) flag")
        T = splitting_involution(group)
        components = indecomposable_components(G, max_nodes)
        require(row["decomposable"] == (T is not None) == (components >= 2) and row["eichler_components"] == components,
                "decomposability")
        if row["decomposable"]:
            inv = matrix(row["decomposition"]["involution"], "involution", 4)
            require(inv in group and multiply(inv, inv) == I4 and inv not in (I4, MINUS_I4)
                    and all((x + y) % 2 == 0 for rT, rI in zip(inv, I4) for x, y in zip(rT, rI)), "stored involution")
            grams = row["decomposition"]["summand_grams"]
            require(prod(det4(g) if len(g) == 4 else det_small(g) for g in grams) == det4(G), "summand determinants")
        else:
            indec += 1
            status = check_decisions(G, group, row["parity"], max_nodes, enumeration_nodes, external_ok=False)
            holds += status == "holds"
            fails += status == "fails"
            unresolved += status == "unresolved"
            if status == "fails":
                failing.append(det4(G))
            for k, v in row["parity"]["first_nonzero_norm4_distribution"].items():
                firsts[int(k)] += v
            for k, v in row["parity"]["zero_levels"].items():
                levels[k] += v
        if progress:
            progress(f"  replayed class {i + 1}/{len(rows)} det {det4(G)} {time.perf_counter() - start:.0f}s")
    s = doc["summary"]
    require(s["classes"] == len(rows) and s["indecomposable_classes"] == indec and s["P_holds"] == holds
            and s["P_fails"] == fails and s["P_unresolved"] == unresolved and sorted(failing) == s["failing_determinants"],
            "wrong decided summary")
    require({int(k): v for k, v in s["first_nonzero_norm4_distribution"].items()} == dict(firsts), "wrong first-nonzero summary")
    require(s["zero_levels"] == dict(levels), "wrong level summary")
    return {"verdict": "verified", **listing, "indecomposable_classes": indec, "P_holds": holds, "P_fails": fails,
            "P_unresolved": unresolved, "failing_determinants": sorted(failing), "seconds": time.perf_counter() - start}


def det_small(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if n == 3:
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    return det4(M)


# ----------------------------------------------------------------------------
# Task A replay
# ----------------------------------------------------------------------------

def replay_task_a(doc, max_nodes, enumeration_nodes, progress=None):
    require(doc["schema"] == "work7-rank4-parity-v1-task-a", "wrong Task A schema")
    start = time.perf_counter()
    census_doc = doc["census"]
    raw = Path(census_doc["source"]).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == census_doc["source_sha256"], "census hash")
    census = json.loads(raw)
    rows = census_doc["lattices"]
    require(len(rows) == len(census["lattices"]) == 169, "census coverage")
    indec = 0
    for row, L in zip(rows, census["lattices"]):
        G, group = L["G"], L["automorphisms"]
        require(row["G"] == G and row["determinant"] == det4(G) == L["determinant"] and row["aut_order"] == len(group)
                and row["even"] == even(G) and row["det_mod_8"] == det4(G) % 8, "census row")
        T = splitting_involution(group)
        components = indecomposable_components(G, max_nodes)
        require(row["decomposable"] == (T is not None) == (components >= 2) and row["eichler_components"] == components, "census decomposability")
        pairs = [p for p in L["pairs"] if any(p["a"]) and any(p["b"]) and dot(p["a"], p["b"]) % 2 == 0]
        zeros = [p for p in pairs if p["verdict"] == "proved_zero"]
        require(len(pairs) == 105 and row["even_zeros"] == len(zeros) and row["P"] == ("fails" if zeros else "holds"), "census P status")
        for p in zeros:
            require(symmetry_data(G, p["a"], p["b"], p["evidence"]["T"])["epsilon"] == 1, "census zero witness")
        if even(G):
            require([(z["a"], z["b"]) for z in row["zeros"]] == [(p["a"], p["b"]) for p in zeros], "census zero list")
            for z in row["zeros"]:
                f2_check(G, z["a"], z["b"], z["f2"])
        indec += not row["decomposable"]
    require(indec == 27 == census_doc["summary"]["indecomposable"], "27 indecomposable census classes expected")
    if progress:
        progress("  census replayed")
    eis = {}
    for key in ("eisenstein_task_b", "eisenstein_task_d"):
        e = doc[key]
        raw = Path(e["source"]).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == e["source_sha256"], "eisenstein hash")
        src = json.loads(raw)
        require(len(e["lattices"]) == len(src["representatives"]) == e["representatives"], "eisenstein coverage")
        for n, (row, r) in enumerate(zip(e["lattices"], src["representatives"])):
            m = src["members"][r["member_index"]]
            G = m["L_reduced"]
            require(row["G"] == G and row["a"] == m["a"] and row["b"] == m["b"] and row["level"] == r["level"]
                    and row["determinant"] == det4(G) == r["primitive_determinant"] and row["even"] == even(G), "eisenstein row")
            group = box_isometries(G, G, max_nodes)
            T = splitting_involution(group)
            components = indecomposable_components(G, max_nodes)
            require(row["aut_order"] == len(group) and row["decomposable"] == (T is not None) == (components >= 2)
                    and row["eichler_components"] == components, "eisenstein decomposability")
            if even(G):
                f2_check(G, m["a"], m["b"], row["f2"])
            if progress and (n + 1) % 200 == 0:
                progress(f"  {key} {n + 1}/{len(e['lattices'])} {time.perf_counter() - start:.0f}s")
        cands = [row for row in e["lattices"] if not row["decomposable"] and not row["even"] and row["det_mod_4"] != 0]
        for lvl in ("1", "2"):
            best = min((row for row in cands if row["level"] == int(lvl)), key=lambda r: r["determinant"], default=None)
            stored = e["smallest_indecomposable_odd_lattice_det_not_div_4_by_level"][lvl]
            require((best is None) == (stored is None) and (best is None or stored["determinant"] == best["determinant"]), "smallest by level")
            eis[f"{key}_smallest_level_{lvl}"] = None if best is None else best["determinant"]
    lv = doc["level4"]
    raw = Path(lv["pairs_source"]).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == lv["pairs_sha256"], "pairs hash")
    pairs = json.loads(raw)["pairs"]
    members = lv["pair_members"]
    require([(m["pair"], m["which"], m["G"]) for m in members] == [(p["name"], w, p[w]) for p in pairs for w in ("G1", "G2")], "pair coverage")
    statuses = Counter()
    for m in members:
        G = m["G"]
        group = box_isometries(G, G, max_nodes)
        require(m["aut_order"] == len(group) and m["decomposable"] == (splitting_involution(group) is not None), "pair member group")
        statuses[check_decisions(G, group, m["parity"], max_nodes, enumeration_nodes)] += 1
        require(m["P"] == m["parity"]["P"], "pair member status")
        if progress:
            progress(f"  pair member {m['pair']} {m['which']} replayed {time.perf_counter() - start:.0f}s")
    l4 = Counter()
    for m in lv["level4_lattices"]:
        G = m["G"]
        group = box_isometries(G, G, max_nodes)
        require(m["aut_order"] == len(group) and m["decomposable"] == (splitting_involution(group) is not None), "level-4 group")
        status = check_decisions(G, group, m["parity"], max_nodes, enumeration_nodes)
        require(status == "fails" == m["P"] or (m["P"] == "fails" and status == "fails"), "level-4 lattice must fail P")
        l4[status] += 1
        if progress:
            progress(f"  level-4 det {m['determinant']} replayed {time.perf_counter() - start:.0f}s")
    return {"verdict": "verified", "census_classes": 169, "census_indecomposable": 27, **eis,
            "pair_members": dict(statuses), "level4_lattices": dict(l4), "seconds": time.perf_counter() - start}


# ----------------------------------------------------------------------------
# (C1) counterexample certificates
# ----------------------------------------------------------------------------

def replay_certificate(cert, max_nodes, enumeration_nodes):
    require(cert["schema"] == "work7-rank4-parity-v1-c1-counterexample", "wrong certificate schema")
    G = matrix(cert["G"], "G", 4)
    ldl(G)
    a, b = vector(cert["a"], "a", 4), vector(cert["b"], "b", 4)
    det = det4(G)
    require(cert["determinant"] == det and det % 2 == 1 and even(G) and cert["even"], "not an even lattice of odd determinant")
    g = 0
    for row in G:
        for x in row:
            g = gcd(g, abs(x))
    require(g == 1, "not primitive")
    require(all(v in (0, 1) for v in a + b) and any(a) and any(b) and dot(a, b) % 2 == 0, "not an even characteristic")
    f2_check(G, a, b, cert["f2"])
    group = box_isometries(G, G, max_nodes)
    require(cert["automorphisms"] == group and cert["aut_order"] == len(group), "full group")
    duals = [transpose(inverse(T)) for T in group]
    indices, signs = stabilizer(G, a, b, group, duals)
    require([s["T"] for s in cert["stabilizer"]] == [group[i] for i in indices] and [s["epsilon"] for s in cert["stabilizer"]] == signs, "stabilizer")
    for s in cert["stabilizer"]:
        T = s["T"]
        P, k = T, 1
        while P != I4:
            P = multiply(P, T)
            k += 1
            require(k <= 64, "order")
        require(s["order"] == k, "stabilizer element order")
    require(splitting_involution(group) is None and indecomposable_components(G, max_nodes) == 1 and cert["indecomposable"], "decomposable")
    cutoff = cert["identity"]["modular_cutoff"]
    require(cutoff["G"] == G and cutoff["a"] == a and cutoff["b"] == b, "cutoff for another characteristic")
    B = modular_check(G, a, b, cutoff["modular"])
    require(cutoff["run"]["norm4_bound"] >= B and cutoff["evidence"]["enumeration_complete"] is True, "cutoff incomplete")
    got = shells(G, a, b, cutoff["run"]["norm4_bound"], enumeration_nodes)
    require(got == cutoff["evidence"]["enumeration"]["shells"] and all(s["signed_coefficient"] == 0 for s in got), "cutoff coefficients")
    kernel = {k: cert[k] for k in ("kernel_basis", "G0", "d0", "complement_t", "classes", "coset_lattices", "level")}
    kernel.update({"aut_L0_order": cert["kernel_automorphisms"]["order"],
                   "elements_mapping_v1_to_v2": cert["kernel_automorphisms"]["elements_mapping_v1_to_v2"],
                   "elements_swapping_v1_v2": cert["kernel_automorphisms"]["elements_swapping_v1_v2"]})
    level = check_kernel(G, a, b, kernel, max_nodes)
    witnesses = [s for s in cert["stabilizer"] if s["epsilon"] == 1]
    require((level == 1) == bool(witnesses), "level and witnesses")
    if witnesses:
        W = cert["identity"]["symmetry_witness"]["T"]
        require(W in group and symmetry_data(G, a, b, W)["epsilon"] == 1, "symmetry witness")
    return {"verdict": "verified_c1_counterexample", "determinant": det, "level": level, "aut_order": len(group),
            "cutoff": B, "witness_orders": sorted({s["order"] for s in witnesses})}


def gcd(x, y):
    while y:
        x, y = y, x % y
    return x


def replay_pair_members(doc, max_nodes, enumeration_nodes):
    """P certificates of the eighteen isospectral-pair members (examples/rank4-parity)."""
    require(doc["schema"] == "work7-rank4-parity-v1-pair-members", "wrong schema")
    raw = Path(doc["pairs_source"]).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == doc["pairs_sha256"], "pairs hash")
    pairs = json.loads(raw)["pairs"]
    members = doc["members"]
    require([(m["pair"], m["which"], m["G"]) for m in members] == [(p["name"], w, p[w]) for p in pairs for w in ("G1", "G2")], "pair coverage")
    statuses = []
    for m in members:
        group = box_isometries(m["G"], m["G"], max_nodes)
        require(m["decomposable"] == (splitting_involution(group) is not None), "decomposability")
        statuses.append({"pair": m["pair"], "which": m["which"], "determinant": det4(m["G"]),
                         "P": check_decisions(m["G"], group, m["parity"], max_nodes, enumeration_nodes)})
        require(statuses[-1]["P"] == m["P"], "stored status")
    return {"verdict": "verified", "members": statuses}


def guarded(fn):
    """Malformed documents are rejected as invalid certificates; node limits still raise."""
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except InvalidCertificate:
            raise
        except (KeyError, TypeError, IndexError, ValueError, ZeroDivisionError, AttributeError, OSError) as error:
            raise InvalidCertificate("malformed document: " + repr(error)) from error
    wrapper.__name__ = fn.__name__
    wrapper.__doc__ = fn.__doc__
    return wrapper


replay_decided = guarded(replay_decided)
replay_task_a = guarded(replay_task_a)
replay_certificate = guarded(replay_certificate)
replay_pair_members = guarded(replay_pair_members)
replay_class_list = guarded(replay_class_list)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("classes", "task-a", "certificates", "pair-members"))
    parser.add_argument("document", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-nodes", type=int, default=5_000_000)
    parser.add_argument("--enumeration-nodes", type=int, default=400_000_000)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--fast-box", action="store_true", help="use the outer-box backend for every bound")
    args = parser.parse_args()
    FAST_BOX[0] = args.fast_box
    import sys
    progress = (lambda msg: print(msg, file=sys.stderr, flush=True)) if args.progress else None
    if args.output:
        require(not args.output.exists(), "choose a fresh output path")
    raw = args.document.read_bytes()
    doc = json.loads(raw)
    if args.kind == "classes":
        result = replay_decided(doc, args.max_nodes, args.enumeration_nodes, progress)
    elif args.kind == "task-a":
        result = replay_task_a(doc, args.max_nodes, args.enumeration_nodes, progress)
    elif args.kind == "pair-members":
        result = replay_pair_members(doc, args.max_nodes, args.enumeration_nodes)
    else:
        result = {"verdict": "verified", "certificates": [replay_certificate(c, args.max_nodes, args.enumeration_nodes)
                                                           for c in doc["certificates"]]}
    result["input"] = str(args.document)
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return 0 if result["verdict"].startswith("verified") else 2


if __name__ == "__main__":
    raise SystemExit(main())
