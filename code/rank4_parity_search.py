"""Rank-four parity property P: tasks A-E of Fable-rank4-parity-conjecture.md.

P(L): every even characteristic (a, b nonzero mod 2, a.b even) of the positive
definite integral rank-four lattice L has a theta function that is not
identically zero.  Decisions are exact: a nonzero is one complete nonzero
signed shell coefficient (exact_theta.complete_shells), a zero is an
epsilon-one stabilizer element (PROOFS.md section 2) or a complete modular
cutoff (RANK4-FEASIBILITY.md section 2); interrupted enumerations are
reported, never counted.  Standard library only.  Report:
RANK4-PARITY-CONJECTURE.md.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import ceil, floor, prod
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells, dot,
                         identity, integer, inverse, ldl, matrix, matvec, multiply,
                         require, symmetry_data, transpose, vector)
from rank4_feasibility import (action_tables, candidates as feasibility_candidates,
                               cutoff_build, determinant, isometries, modular_data)
from rank4_level4_search import (class_index, class_vector, gram_of, hnf_basis,
                                 kernel_basis, lll_reduce, orbits_on_classes,
                                 overlattice_basis, reduced_lattice,
                                 transport_characteristic)
from lattice_census import sphere


def _source(path):
    """Resolve a hash-bound research-workspace path (examples/, results/) in the release layout (data/, reports/)."""
    p = Path(path)
    if p.exists():
        return p
    s = str(path).replace("\\", "/")
    for old, new in (("examples/", "data/"), ("results/", "reports/")):
        if s.startswith(old):
            return Path(__file__).resolve().parent.parent / (new + s[len(old):])
    return p

SCHEMA = "work7-rank4-parity-v1"
RATIO = Fraction(3, 4)
BITS = [list(v) for v in product((0, 1), repeat=4)]
EVEN_CHARACTERISTICS = [(a, b) for a in BITS if any(a)
                        for b in BITS if any(b) and dot(a, b) % 2 == 0]   # 105 pairs
I4 = identity(4)
MINUS_I4 = [[-x for x in row] for row in I4]
LADDER = (100, 400, 1600)


def is_even(G):
    return all(G[i][i] % 2 == 0 for i in range(4))


def sha256_json(obj):
    return hashlib.sha256(json.dumps(obj, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


# ----------------------------------------------------------------------------
# decomposability: two criteria
# ----------------------------------------------------------------------------

def splitting_involution(group):
    """An involution T != +-I of Aut(L) with (I+T)/2 integral, or None.

    Lemma.  L is decomposable iff such a T exists.  (=>) If L = A (+) B, the
    map T = (-1 on A, +1 on B) is an isometry of L, an involution, and
    (I+T)/2 is the projection onto B, integral on L.  (<=) P = (I+T)/2 is an
    integral idempotent; T is self-adjoint for the form (T = T^{-1}), so P is
    an orthogonal projection and L = PL (+) (I-P)L with both summands nonzero
    because T != +-I.  Both directions are checked on data by
    decomposition_from_involution and eichler_components.
    """
    for T in group:
        if T == I4 or T == MINUS_I4 or multiply(T, T) != I4:
            continue
        if all((x + y) % 2 == 0 for rT, rI in zip(T, I4) for x, y in zip(rT, rI)):
            return T
    return None


def span_basis(cols):
    """Basis (4 x r matrix of columns) of the Z-span of integer vectors in Z^4."""
    cols = [list(c) for c in cols]
    basis = []
    for i in range(4):
        pivots = [c for c in cols if c[i] != 0]
        rest = [c for c in cols if c[i] == 0]
        while len(pivots) > 1:
            pivots.sort(key=lambda c: abs(c[i]))
            p = pivots[0]
            new = [p]
            for c in pivots[1:]:
                q = c[i] // p[i]
                c = [x - q * y for x, y in zip(c, p)]
                (new if c[i] != 0 else rest).append(c)
            pivots = new
        if pivots:
            p = pivots[0]
            basis.append([-x for x in p] if p[i] < 0 else p)
        cols = rest
    return transpose(basis)


def decomposition_from_involution(G, T):
    """Explicit orthogonal decomposition L = PL (+) (I-P)L from a splitting involution."""
    P = [[(x + y) // 2 for x, y in zip(rT, rI)] for rT, rI in zip(T, I4)]
    Q = [[y - x for x, y in zip(rP, rI)] for rP, rI in zip(P, I4)]
    summands = []
    for M in (P, Q):
        basis = span_basis([c for c in transpose(M) if any(c)])
        summands.append({"basis": basis, "gram": gram_of(G, basis), "rank": len(basis[0])})
    B = [r1 + r2 for r1, r2 in zip(summands[0]["basis"], summands[1]["basis"])]
    require(abs(determinant(B)) == 1, "summand bases do not form a basis of L")
    full = gram_of(G, B)
    r = summands[0]["rank"]
    require(all(full[i][j] == 0 for i in range(r) for j in range(r, 4)), "summands not orthogonal")
    return {"involution": T, "summands": summands, "basis_change": B}


def eichler_components(G, max_nodes=5_000_000):
    """Components of the graph of indecomposable vectors of norm <= max G_ii.

    A vector v is decomposable if v = u + w with u, w nonzero and orthogonal,
    equivalently some u has 0 < u.u < v.v and u.v = u.u.  Vectors from
    different components are pairwise orthogonal, so the sublattices they
    generate are orthogonal; every basis vector is a sum of indecomposable
    vectors of no larger norm, so these sublattices sum to L.  Hence the
    number of components is at least two iff L is decomposable (if L = A (+) B,
    every indecomposable vector lies in A or in B, and the two sides are
    orthogonal, so the graph is disconnected).
    """
    R = max(G[i][i] for i in range(4))
    vecs = [v for v in sphere(G, R, max_nodes) if any(v)]
    norms = {v: dot(v, matvec(G, v)) for v in vecs}
    by_norm = defaultdict(list)
    for v in vecs:
        by_norm[norms[v]].append(v)
    indecomposable = []
    for v in vecs:
        if v < tuple(-x for x in v):
            continue                       # one of each antipodal pair
        Gv = matvec(G, v)
        decomposable = False
        for n in sorted(by_norm):
            if n >= norms[v]:
                break
            if any(dot(u, Gv) == n for u in by_norm[n]):
                decomposable = True
                break
        if not decomposable:
            indecomposable.append(v)
    parent = list(range(len(indecomposable)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    Gi = [matvec(G, v) for v in indecomposable]
    for i in range(len(indecomposable)):
        for j in range(i + 1, len(indecomposable)):
            if dot(indecomposable[i], Gi[j]) != 0:
                parent[find(i)] = find(j)
    components = len({find(i) for i in range(len(indecomposable))})
    return {"norm_bound": R, "indecomposable_vectors_up_to_sign": len(indecomposable),
            "components": components}


def decomposability(G, group, max_nodes=5_000_000):
    """Both criteria, required to agree."""
    T = splitting_involution(group)
    eich = eichler_components(G, max_nodes)
    require((T is not None) == (eich["components"] >= 2), "decomposability criteria disagree")
    out = {"decomposable": T is not None, "eichler": eich}
    if T is not None:
        out["decomposition"] = decomposition_from_involution(G, T)
    return out


# ----------------------------------------------------------------------------
# the quadratic space L/2L of an even lattice
# ----------------------------------------------------------------------------

def mod2_rank_and_kernel(M):
    """Rank over F2 and a basis of the kernel of the 4 x 4 integer matrix M mod 2."""
    rows = [[x % 2 for x in row] for row in M]
    pivcols, r = [], 0
    for c in range(4):
        piv = next((i for i in range(r, 4) if rows[i][c]), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(4):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        pivcols.append(c)
        r += 1
    kernel = []
    for f in [c for c in range(4) if c not in pivcols]:
        v = [0] * 4
        v[f] = 1
        for i, c in enumerate(pivcols):
            v[c] = rows[i][f]
        kernel.append(v)
    return r, kernel


def solve_mod2(M, b):
    """All x in F2^4 with M x = b mod 2."""
    return [list(x) for x in BITS if all((y - z) % 2 == 0 for y, z in zip(matvec(M, x), b))]


def quadratic_value(G, x):
    """Q(x) = x^t G x / 2 mod 2 for an even G."""
    n = dot(x, matvec(G, x))
    require(n % 2 == 0, "Q needs an even lattice")
    return (n // 2) % 2


def f2_data(G, a, b):
    """(a, c) geometry in (L/2L, Q) for an even lattice, radical included."""
    require(is_even(G), "F2 geometry is defined for even lattices")
    rank, radical = mod2_rank_and_kernel(G)
    rad_span = sorted({tuple(sum(k * r[i] for k, r in zip(coeffs, radical)) % 2 for i in range(4))
                       for coeffs in product((0, 1), repeat=len(radical))})
    singular = [x for x in BITS if any(x) and quadratic_value(G, x) == 0]
    cs = solve_mod2(G, b)
    Q_rad = [quadratic_value(G, list(r)) for r in rad_span]
    out = {"rank_mod_2": rank, "radical_dimension": 4 - rank, "radical": [list(r) for r in rad_span],
           "Q_on_radical": Q_rad, "radical_totally_singular": all(q == 0 for q in Q_rad),
           "nonzero_singular_vectors": len(singular),
           "type": ("plus" if len(singular) == 9 else "minus" if len(singular) == 5 else None) if rank == 4 else "degenerate",
           "a_in_radical": list(a) in [list(r) for r in rad_span],
           "Q_a": quadratic_value(G, a), "b_in_image_mod_2": bool(cs), "c_solutions": cs,
           "pairing_a_b_mod_2": dot(a, b) % 2}
    if cs:
        out["Q_c_values"] = sorted({quadratic_value(G, c) for c in cs})
        out["a_equals_some_c"] = any(list(a) == c for c in cs)
        out["span_a_c_totally_singular"] = (out["Q_a"] == 0 and out["Q_c_values"] == [0]
                                            and not out["a_equals_some_c"])
        out["a_in_radical_plus_c"] = any(all((x - y - r) % 2 == 0 for x, y, r in zip(a, c, rad))
                                         for c in cs for rad in rad_span)
    return out


def arf_type(G):
    if not is_even(G):
        return None
    rank, _ = mod2_rank_and_kernel(G)
    singular = sum(1 for x in BITS if any(x) and quadratic_value(G, x) == 0)
    if rank < 4:
        return {"type": "degenerate", "rank_mod_2": rank, "nonzero_singular_vectors": singular}
    return {"type": "plus" if singular == 9 else "minus", "rank_mod_2": 4, "nonzero_singular_vectors": singular}


# ----------------------------------------------------------------------------
# deciding P for one lattice
# ----------------------------------------------------------------------------

def matrix_order(T, limit=64):
    P = T
    for k in range(1, limit + 1):
        if P == I4:
            return k
        P = multiply(P, T)
    raise InvalidCertificate("order exceeds the limit")


def stabilizer_signs(G, a, b, group, tables=None):
    """Indices of the stabilizing elements of (a, b) in group and their epsilon values."""
    if tables is None:
        tables = action_tables(group)
    ai, bi = BITS.index(list(a)), BITS.index(list(b))
    indices, signs = [], []
    for i, (primal, dual) in enumerate(tables):
        if primal[ai] == list(a) and all((v - w) % 2 == 0 for v, w in zip(dual[bi], b)):
            indices.append(i)
            signs.append(dot(a, [(v - w) // 2 for v, w in zip(dual[bi], b)]) % 2)
    return indices, signs


def first_shell(shells):
    s = next((s for s in shells if s["vector_count"]), None)
    return None if s is None else {"norm4": s["norm4"], "vector_count": s["vector_count"],
                                   "signed_coefficient": s["signed_coefficient"]}


def kernel_configuration(G, a, b, max_nodes=2_000_000):
    """Sign kernel L0, the two coset classes and the level of an identity zero."""
    K = kernel_basis(b)
    G0 = gram_of(G, K)
    Kinv = inverse(K)
    j = next(i for i in range(4) if b[i])
    t = [int(i == j) for i in range(4)]
    w1 = matvec(Kinv, a)
    w2 = matvec(Kinv, [p + 2 * q for p, q in zip(a, t)])
    require(all(x.denominator == 1 for x in w1 + w2), "a does not lie in the kernel")
    v1 = [int(x) % 2 for x in w1]
    v2 = [int(x) % 2 for x in w2]
    group0 = isometries(G0, G0, max_nodes)
    perms = [[class_index(matvec(T, class_vector(k))) for k in range(16)] for T in group0]
    k1, k2 = class_index(v1), class_index(v2)
    swapping = [i for i, p in enumerate(perms) if p[k1] == k2 and p[k2] == k1]
    mapping = [i for i, p in enumerate(perms) if p[k1] == k2]
    cosets = {}
    for name, v in (("O_v1", v1), ("O_v2", v2)):
        C = overlattice_basis(v)
        cosets[name] = {"doubled_basis_in_L0": C, "gram_times_4": gram_of(G0, C),
                        "reduced": reduced_lattice(G0, C)}
    out = {"kernel_basis": K, "G0": G0, "d0": determinant(G0), "complement_t": t,
           "classes": {"v1": v1, "v2": v2}, "coset_lattices": cosets,
           "aut_L0_order": len(group0), "orbits_on_classes": orbits_on_classes(group0)[0],
           "elements_mapping_v1_to_v2": len(mapping), "elements_swapping_v1_v2": len(swapping)}
    if swapping:
        out["level"], out["isometry_witness"] = 1, {"kind": "swap in Aut(L0)", "T_on_L0": group0[swapping[0]]}
    elif mapping:
        out["level"], out["isometry_witness"] = 2, {"kind": "map in Aut(L0)", "T_on_L0": group0[mapping[0]]}
    else:
        iso = isometries(cosets["O_v1"]["reduced"]["gram"], cosets["O_v2"]["reduced"]["gram"], max_nodes)
        out["level"] = 3 if iso else 4
        out["isometry_witness"] = {"kind": "isometry of reduced coset lattices", "U": iso[0]} if iso else None
    return out


def decide_characteristic(G, a, b, group, tables, ladder=LADDER, max_nodes=2_000_000,
                          cutoff_vector_limit=5_000_000, cutoff_nodes=200_000_000, want_level=True):
    """Exact decision for one even characteristic; see section 3 of the brief."""
    a, b = list(a), list(b)
    indices, signs = stabilizer_signs(G, a, b, group, tables)
    rec = {"a": a, "b": b, "stabilizer": {"indices": indices, "epsilon": signs}}
    if 1 in signs:
        T = group[indices[signs.index(1)]]
        rec.update(verdict="zero_by_symmetry", witness=T, witness_order=matrix_order(T), level=1)
        return rec
    for bound in ladder:
        try:
            enum = complete_shells(G, a, b, bound, max_nodes)
        except EnumerationLimit as error:
            rec.update(verdict="unresolved", reason=str(error), enumeration_complete=False, norm4_bound=bound)
            return rec
        nz = next((s for s in enum["shells"] if s["signed_coefficient"]), None)
        if nz:
            rec.update(verdict="nonzero", norm4_bound=bound, shells=enum["shells"],
                       first_nonzero={"norm4": nz["norm4"], "coefficient": nz["signed_coefficient"]},
                       minimal_shell=first_shell(enum["shells"]))
            return rec
    item = modular_data(G, a, b)
    B = item["cutoff"]
    det = determinant(G)
    estimate = (3.141592653589793 ** 2 / 2) * B * B / (16 * det ** 0.5)
    rec.update(cutoff=B, estimated_cutoff_vectors=estimate, ladder_exhausted=ladder[-1])
    if estimate > cutoff_vector_limit:
        rec.update(verdict="unresolved", reason=f"all coefficients through norm4 {ladder[-1]} vanish; "
                   f"cutoff {B} estimated at {estimate:.3g} vectors exceeds the limit {cutoff_vector_limit}")
        return rec
    doc = cutoff_build(G, a, b, max_nodes=cutoff_nodes)
    ev = doc["evidence"]
    if not ev["enumeration_complete"]:
        rec.update(verdict="unresolved", reason=ev["reason"], enumeration_complete=False)
        return rec
    if ev["verdict"] == "proved_nonzero":
        rec.update(verdict="nonzero", norm4_bound=B, shells=ev["enumeration"]["shells"],
                   first_nonzero={"norm4": ev["nonzero_norm4"], "coefficient": ev["coefficient"]},
                   minimal_shell=first_shell(ev["enumeration"]["shells"]), cutoff_certificate=doc)
        return rec
    require(ev["verdict"] == "proved_zero_by_modular_cutoff", "unexpected cutoff verdict")
    rec.update(verdict="zero_by_cutoff", cutoff_certificate=doc, minimal_shell=first_shell(ev["enumeration"]["shells"]))
    if want_level:
        rec["kernel"] = kernel_configuration(G, a, b, max_nodes)
        rec["level"] = rec["kernel"]["level"]
    return rec


def decide_parity(G, group=None, max_nodes=2_000_000, ladder=LADDER, cutoff_vector_limit=5_000_000,
                  cutoff_nodes=200_000_000, with_f2=True, want_level=True):
    """All 105 even characteristics of L; P holds iff all are proved nonzero."""
    matrix(G, "G", 4)
    start = time.perf_counter()
    if group is None:
        group = isometries(G, G, max_nodes)
    tables = action_tables(group)
    decisions = [decide_characteristic(G, a, b, group, tables, ladder, max_nodes, cutoff_vector_limit,
                                       cutoff_nodes, want_level) for a, b in EVEN_CHARACTERISTICS]
    verdicts = Counter(d["verdict"] for d in decisions)
    zeros = [d for d in decisions if d["verdict"].startswith("zero")]
    status = "holds" if verdicts["nonzero"] == 105 else "fails" if zeros else "unresolved"
    nonzero = [d for d in decisions if d["verdict"] == "nonzero"]
    even = is_even(G)
    if with_f2 and even:
        for d in decisions:
            d["f2"] = f2_data(G, d["a"], d["b"])
    out = {"G": G, "determinant": determinant(G), "even": even, "aut_order": len(group),
           "P": status, "verdicts": dict(verdicts), "even_zeros": len(zeros),
           "zero_levels": dict(Counter(d.get("level") for d in zeros)),
           "first_nonzero_norm4_distribution": dict(sorted(Counter(d["first_nonzero"]["norm4"] for d in nonzero).items())),
           "largest_first_nonzero_norm4": max((d["first_nonzero"]["norm4"] for d in nonzero), default=None),
           "minimal_shell_decides": sum(1 for d in nonzero if d["first_nonzero"]["norm4"] == d["minimal_shell"]["norm4"]),
           "arf": arf_type(G) if even else None, "decisions": decisions,
           "seconds": time.perf_counter() - start}
    return out


# ----------------------------------------------------------------------------
# certificate of a counterexample to (C1): section 5 of the brief
# ----------------------------------------------------------------------------

def primitive_scaling(G):
    g = 0
    for row in G:
        for x in row:
            g = gcd_int(g, abs(x))
    return g


def gcd_int(a, b):
    while b:
        a, b = b, a % b
    return a


def certify_counterexample(G, a, b, group=None, max_nodes=2_000_000, cutoff_nodes=400_000_000,
                           cutoff_vector_limit=30_000_000):
    """Reduced primitive Gram, characteristic, kernel, cosets, identity proofs, groups, level."""
    start = time.perf_counter()
    matrix(G, "G", 4)
    require(is_even(G) and determinant(G) % 2 == 1, "a (C1) counterexample is even with odd determinant")
    require(primitive_scaling(G) == 1, "the Gram matrix must be primitive")
    U, Gr = lll_reduce(G)
    a, b = transport_characteristic(U, list(a), list(b))
    G = Gr
    det = determinant(G)
    if group is None or U != I4:
        group = isometries(G, G, max_nodes)     # the group of the reduced Gram, recomputed
    tables = action_tables(group)
    indices, signs = stabilizer_signs(G, a, b, group, tables)
    stabilizer = [{"T": group[i], "epsilon": e, "order": matrix_order(group[i])} for i, e in zip(indices, signs)]
    conf = kernel_configuration(G, a, b, max_nodes)
    item = modular_data(G, a, b)
    B = item["cutoff"]
    estimate = (3.141592653589793 ** 2 / 2) * B * B / (16 * det ** 0.5)
    require(estimate <= cutoff_vector_limit, "cutoff too large for this certificate")
    doc = cutoff_build(G, a, b, max_nodes=cutoff_nodes)
    require(doc["evidence"]["enumeration_complete"] and doc["evidence"]["verdict"] == "proved_zero_by_modular_cutoff",
            "the cutoff does not prove the zero")
    decomposition = decomposability(G, group, max_nodes)
    require(not decomposition["decomposable"], "L is decomposable")
    witnesses = [s for s in stabilizer if s["epsilon"] == 1]
    require((conf["level"] == 1) == bool(witnesses), "level and stabilizer disagree")
    identity_proof = {"modular_cutoff": doc, "isometry_witness": conf["isometry_witness"]}
    if witnesses:
        identity_proof["symmetry_witness"] = witnesses[0]
    return {"schema": SCHEMA + "-c1-counterexample", "G": G, "determinant": det, "even": True,
            "lll_transform_from_input": U, "a": a, "b": b, "pairing": dot(a, b),
            "f2": f2_data(G, a, b), "arf": arf_type(G),
            "kernel_basis": conf["kernel_basis"], "G0": conf["G0"], "d0": conf["d0"],
            "complement_t": conf["complement_t"], "classes": conf["classes"],
            "coset_lattices": conf["coset_lattices"], "identity": identity_proof,
            "automorphisms": group, "aut_order": len(group), "stabilizer": stabilizer,
            "element_orders_in_aut": dict(sorted(Counter(matrix_order(T) for T in group).items())),
            "indecomposable": True, "decomposability": decomposition,
            "kernel_automorphisms": {"order": conf["aut_L0_order"], "orbits_on_classes": conf["orbits_on_classes"],
                                     "elements_mapping_v1_to_v2": conf["elements_mapping_v1_to_v2"],
                                     "elements_swapping_v1_v2": conf["elements_swapping_v1_v2"]},
            "level": conf["level"], "seconds": time.perf_counter() - start}


# ----------------------------------------------------------------------------
# Task B generator: the reduced Gram box of RANK4-FEASIBILITY.md section 1, restated
# ----------------------------------------------------------------------------

COVERAGE = ("Every positive definite integral rank-four lattice has a basis v1..v4 obtained by taking a "
            "shortest primitive vector, projecting, repeating, and size-reducing so that all Gram-Schmidt "
            "coefficients satisfy |mu_ij| <= 1/2 (RANK4-FEASIBILITY.md section 1).  Then d1 = G11 = A is a "
            "positive integer, d_(i+1) >= (3/4) d_i, det = d1 d2 d3 d4 >= (3/4)^6 A^4, and G_ii >= A.  The "
            "generator enumerates A with (3/4)^6 A^4 <= D0; for each new basis vector the cross entries G_jk "
            "lie in the integral interval of radius d_j/2 around sum_(l<j) mu_jl mu_kl d_l; the diagonal G_kk "
            "satisfies G_kk >= max(A, q + (3/4) d_(k-1)) with q = sum_j mu_kj^2 d_j and prefix * d_k^r * "
            "(3/4)^(r(r-1)/2) <= D0 for the r remaining pivots; the fourth diagonal solves prefix * (G44 - q) = "
            "det for each admitted det.  The box therefore contains at least one basis of every lattice with "
            "det <= D0.  Restricting to even diagonals discards only bases of odd lattices (an even lattice has "
            "even diagonals in every basis), and the determinant filter discards only bases of lattices outside "
            "the target universe.  Candidates are bases, not classes; exact isometry deduplication follows.")

DET_FILTERS = {"all": lambda d: True, "odd": lambda d: d % 2 == 1,
               "divisible_by_4": lambda d: d % 4 == 0, "not_divisible_by_4": lambda d: d % 4 != 0}


def candidate_grams(D0, even_only=True, det_filter="all", dmin=1, max_nodes=200_000_000):
    """All admitted reduced-basis Gram matrices with dmin <= det <= D0 passing the filters."""
    integer(D0, "D0")
    integer(dmin, "dmin")
    require(1 <= dmin <= D0, "invalid determinant range")
    det_ok = DET_FILTERS[det_filter]
    admitted = [d for d in range(dmin, D0 + 1) if det_ok(d)]
    step = 2 if even_only else 1
    found, nodes = [], 0

    def tick():
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise EnumerationLimit("candidate generation interrupted; partial domain discarded")

    def extend(G):
        size = len(G)
        mu, pivots = ldl(G)
        prefix = prod(pivots)
        remaining = 4 - size

        def crosses(j, entries, coeffs):
            if j == size:
                q = sum(c * c * d for c, d in zip(coeffs, pivots))
                lower = max(G[0][0], ceil(q + RATIO * pivots[-1]))
                if even_only and lower % 2:
                    lower += 1
                if size == 3:
                    diagonals = []
                    for delta in admitted:
                        value = Fraction(delta) / prefix + q
                        if value.denominator == 1 and value >= lower and (not even_only or value % 2 == 0):
                            diagonals.append(int(value))
                else:
                    diagonals = []
                    diagonal = lower
                    while prefix * (Fraction(diagonal) - q) ** remaining * RATIO ** (remaining * (remaining - 1) // 2) <= D0:
                        tick()
                        diagonals.append(diagonal)
                        diagonal += step
                for diagonal in diagonals:
                    tick()
                    row = [list(old) + [entry] for old, entry in zip(G, entries)]
                    row.append(entries + [diagonal])
                    if size == 3:
                        found.append(row)
                    else:
                        extend(row)
                return
            center = sum(mu[j][k] * coeffs[k] * pivots[k] for k in range(j))
            for entry in range(ceil(center - pivots[j] / 2), floor(center + pivots[j] / 2) + 1):
                tick()
                crosses(j + 1, entries + [entry], coeffs + [(entry - center) / pivots[j]])

        crosses(0, [], [])

    A = step
    while RATIO ** 6 * A ** 4 <= D0:
        tick()
        extend([[A]])
        A += step
    return sorted(found)


def theta_invariant(G, norm_bound=12, max_nodes=2_000_000):
    """Vector counts of the lattice by norm <= norm_bound: an isometry invariant."""
    enum = complete_shells(G, [0] * 4, [0] * 4, 4 * norm_bound, max_nodes)
    return [[s["norm4"] // 4, s["vector_count"]] for s in enum["shells"]]


def class_list(D0, even_only=True, det_filter="odd", dmin=1, max_nodes=5_000_000, progress=None,
               odd_lattices_only=False):
    """Isometry classes in the universe, with generation transports for replay."""
    start = time.perf_counter()
    forms = candidate_grams(D0, even_only, det_filter, dmin)
    if odd_lattices_only:
        forms = [G for G in forms if not is_even(G)]
    t_gen = time.perf_counter() - start
    reps, generation = [], []
    by_key = defaultdict(list)
    searches = 0
    for n, G in enumerate(forms):
        det = determinant(G)
        inv = theta_invariant(G)
        key = (det, json.dumps(inv))
        match = None
        for i in by_key[key]:
            searches += 1
            maps = isometries(G, reps[i]["G"], max_nodes)
            if maps:
                match = (i, maps[0])
                break
        if match is None:
            match = (len(reps), I4)
            reps.append({"G": G, "determinant": det, "invariant": inv})
            by_key[key].append(len(reps) - 1)
        generation.append({"candidate": n, "representative": match[0], "T": match[1]})
        if progress and (n + 1) % 1000 == 0:
            progress(f"  candidates {n + 1}/{len(forms)} classes {len(reps)} {time.perf_counter() - start:.0f}s")
    universe = (f"all isometry classes of positive definite integral rank-four lattices with "
                f"{dmin} <= det <= {D0}, det {det_filter}, "
                + ("even lattices only" if even_only else "odd lattices only" if odd_lattices_only else "all parities"))
    return {"schema": SCHEMA + "-class-list", "universe": universe, "determinant_bound": D0, "determinant_minimum": dmin,
            "even_only": even_only, "odd_lattices_only": odd_lattices_only, "determinant_filter": det_filter,
            "coverage_argument": COVERAGE, "candidate_bases": forms, "candidates": len(forms),
            "classes": len(reps), "representatives": reps, "generation": generation,
            "invariant": "determinant and the vector counts of all norms <= 12 (complete LDL enumeration)",
            "dedup": "exact isometry: complete basis-image search (rank4_feasibility.isometries) between candidates "
                     "with equal invariants; a differing invariant is an exact proof of non-isometry",
            "isometry_searches": searches,
            "seconds": {"generation": t_gen, "total": time.perf_counter() - start}}


def generator_checks(census_path="examples/census-rank4/det24-bound32.json", max_nodes=5_000_000):
    """The fresh generator against the accepted census and the frozen det-24 candidate domain."""
    start = time.perf_counter()
    unrestricted = candidate_grams(24, even_only=False, det_filter="all")
    require(unrestricted == feasibility_candidates(24), "unrestricted generator differs from rank4_feasibility.candidates")
    raw = _source(census_path).read_bytes()
    census = json.loads(raw)
    target = [(L["G"], L["determinant"]) for L in census["lattices"] if is_even(L["G"]) and L["determinant"] % 2 == 1]
    mine = class_list(24, True, "odd", 1, max_nodes)
    matched = []
    for G, det in target:
        hits = [i for i, r in enumerate(mine["representatives"]) if r["determinant"] == det and isometries(G, r["G"], max_nodes)]
        require(len(hits) == 1, "census class without a unique generated class")
        matched.append(hits[0])
    require(sorted(matched) == list(range(mine["classes"])), "generated class not in the census")
    return {"unrestricted_det24_candidates": len(unrestricted), "matches_rank4_feasibility_candidates": True,
            "census_sha256": hashlib.sha256(raw).hexdigest(),
            "census_even_odd_det_classes": [{"determinant": d, "G": G} for G, d in target],
            "generated_even_odd_det_classes": mine["classes"], "generated_candidates": mine["candidates"],
            "bijection_with_census": True, "determinants": dict(sorted(Counter(d for _, d in target).items())),
            "seconds": time.perf_counter() - start}


# ----------------------------------------------------------------------------
# Tasks B, C, E: class lists with P decisions
# ----------------------------------------------------------------------------

def decide_class_list(classes, max_nodes=5_000_000, cutoff_vector_limit=5_000_000, progress=None,
                      decide_decomposable=False):
    start = time.perf_counter()
    rows = []
    for i, rep in enumerate(classes["representatives"]):
        G = rep["G"]
        group = isometries(G, G, max_nodes)
        dec = decomposability(G, group, max_nodes)
        det = rep["determinant"]
        row = {"index": i, "G": G, "determinant": det, "det_mod_4": det % 4, "det_mod_8": det % 8,
               "even": is_even(G), "aut_order": len(group), "automorphisms": group,
               "element_orders": dict(sorted(Counter(matrix_order(T) for T in group).items())),
               "has_square_root_of_minus_identity": any(multiply(T, T) == MINUS_I4 for T in group),
               "decomposable": dec["decomposable"], "eichler_components": dec["eichler"]["components"],
               "arf": arf_type(G)}
        if dec["decomposable"]:
            row["decomposition"] = {"involution": dec["decomposition"]["involution"],
                                    "summand_grams": [s["gram"] for s in dec["decomposition"]["summands"]]}
        if not dec["decomposable"] or decide_decomposable:
            row["parity"] = decide_parity(G, group, max_nodes, LADDER, cutoff_vector_limit)
        rows.append(row)
        if progress:
            progress(f"  class {i + 1}/{classes['classes']} det {det} aut {len(group)} "
                     f"{'decomposable' if dec['decomposable'] else row['parity']['P']} {time.perf_counter() - start:.0f}s")
    indec = [r for r in rows if not r["decomposable"]]
    summary = {"classes": len(rows), "decomposable_classes": len(rows) - len(indec), "indecomposable_classes": len(indec),
               "P_holds": sum(1 for r in indec if r["parity"]["P"] == "holds"),
               "P_fails": sum(1 for r in indec if r["parity"]["P"] == "fails"),
               "P_unresolved": sum(1 for r in indec if r["parity"]["P"] == "unresolved"),
               "failing_determinants": sorted(r["determinant"] for r in indec if r["parity"]["P"] == "fails"),
               "first_nonzero_norm4_distribution": dict(sorted(Counter(
                   int(k) for r in indec for k, v in r["parity"]["first_nonzero_norm4_distribution"].items() for _ in range(v)).items())),
               "largest_first_nonzero_norm4": max((r["parity"]["largest_first_nonzero_norm4"] or 0 for r in indec), default=None),
               "zero_levels": dict(Counter(lvl for r in indec for lvl, n in r["parity"]["zero_levels"].items() for _ in range(n))),
               "by_type": {t: {"classes": sum(1 for r in indec if (r["arf"] or {}).get("type") == t),
                               "P_holds": sum(1 for r in indec if (r["arf"] or {}).get("type") == t and r["parity"]["P"] == "holds"),
                               "P_fails": sum(1 for r in indec if (r["arf"] or {}).get("type") == t and r["parity"]["P"] == "fails")}
                           for t in ("plus", "minus", "degenerate")},
               "seconds": time.perf_counter() - start}
    return {"schema": SCHEMA + "-decided-classes", "classes_doc": classes, "lattices": rows, "summary": summary}


def task_b(D0, max_nodes=5_000_000, cutoff_vector_limit=5_000_000, progress=None):
    classes = class_list(D0, True, "odd", 1, max_nodes, progress)
    if progress:
        progress(f"class list: {classes['candidates']} candidates, {classes['classes']} classes, "
                 f"{classes['seconds']['total']:.1f}s")
    doc = decide_class_list(classes, max_nodes, cutoff_vector_limit, progress)
    doc["task"] = "B"
    doc["conjecture_C1"] = ("refuted" if doc["summary"]["P_fails"] else
                            "holds through D0" if doc["summary"]["P_unresolved"] == 0 else "unresolved classes remain")
    return doc


def task_c(D0=64, max_nodes=5_000_000, cutoff_vector_limit=5_000_000, progress=None):
    classes = class_list(D0, True, "divisible_by_4", 1, max_nodes, progress)
    doc = decide_class_list(classes, max_nodes, cutoff_vector_limit, progress)
    doc["task"] = "C"
    return doc


def task_e_estimate(bounds=(24, 48, 72, 96), max_nodes=5_000_000):
    """Measured generation sizes for all-parity candidates, to extrapolate the odd-lattice listing."""
    rows = []
    for H in bounds:
        t0 = time.perf_counter()
        forms = candidate_grams(H, even_only=False, det_filter="not_divisible_by_4")
        odd = [G for G in forms if not is_even(G)]
        rows.append({"determinant_bound": H, "candidates": len(forms), "odd_lattice_candidates": len(odd),
                     "generation_seconds": time.perf_counter() - t0})
    return rows


def task_e(dmin=25, D0=168, max_nodes=5_000_000, cutoff_vector_limit=5_000_000, progress=None):
    classes = class_list(D0, False, "not_divisible_by_4", dmin, max_nodes, progress, odd_lattices_only=True)
    if progress:
        progress(f"class list: {classes['candidates']} candidates, {classes['classes']} classes, "
                 f"{classes['seconds']['total']:.1f}s")
    doc = decide_class_list(classes, max_nodes, cutoff_vector_limit, progress)
    doc["task"] = "E"
    return doc


# ----------------------------------------------------------------------------
# Task A: the evidence table
# ----------------------------------------------------------------------------

def task_a_census(path="examples/census-rank4/det24-bound32.json", max_nodes=5_000_000):
    raw = _source(path).read_bytes()
    doc = json.loads(raw)
    require(doc["schema"] == "work7-rank4-census-v1" and doc["summary"]["lattice_classes"] == 169, "unexpected census")
    rows = []
    for idx, L in enumerate(doc["lattices"]):
        G, det, group = L["G"], L["determinant"], L["automorphisms"]
        require(determinant(G) == det, "census determinant")
        dec = decomposability(G, group, max_nodes)
        pairs = [p for p in L["pairs"] if any(p["a"]) and any(p["b"]) and dot(p["a"], p["b"]) % 2 == 0]
        require(len(pairs) == 105 and all(p["verdict"] in ("proved_zero", "proved_nonzero") for p in pairs),
                "census even characteristics are not all decided")
        zeros = [p for p in pairs if p["verdict"] == "proved_zero"]
        nonzero = [p for p in pairs if p["verdict"] == "proved_nonzero"]
        row = {"index": idx, "G": G, "determinant": det, "det_mod_4": det % 4, "det_mod_8": det % 8,
               "even": is_even(G), "aut_order": len(group),
               "decomposable": dec["decomposable"], "eichler_components": dec["eichler"]["components"],
               "has_square_root_of_minus_identity": any(multiply(T, T) == MINUS_I4 for T in group),
               "P": "fails" if zeros else "holds", "even_zeros": len(zeros),
               "zero_witness_orders": dict(sorted(Counter(matrix_order(p["evidence"]["T"]) for p in zeros).items())),
               "first_nonzero_norm4_distribution": dict(sorted(Counter(p["evidence"]["norm4"] for p in nonzero).items()))}
        if dec["decomposable"]:
            row["summand_grams"] = [s["gram"] for s in dec["decomposition"]["summands"]]
        if is_even(G):
            row["arf"] = arf_type(G)
            row["zeros"] = [{"a": p["a"], "b": p["b"], "witness_order": matrix_order(p["evidence"]["T"]),
                             "f2": f2_data(G, p["a"], p["b"])} for p in zeros]
            row["totally_singular_pairs"] = sum(1 for p in pairs if f2_data(G, p["a"], p["b"]).get("span_a_c_totally_singular"))
        rows.append(row)
    indec = [r for r in rows if not r["decomposable"]]
    summary = {"classes": len(rows), "indecomposable": len(indec),
               "every_decomposable_class_fails_P": all(r["P"] == "fails" for r in rows if r["decomposable"]),
               "indecomposable_holding": dict(Counter(f"{'even' if r['even'] else 'odd'},det%4={r['det_mod_4']}" for r in indec if r["P"] == "holds")),
               "indecomposable_failing": dict(Counter(f"{'even' if r['even'] else 'odd'},det%4={r['det_mod_4']}" for r in indec if r["P"] == "fails")),
               "even_classes_det_mod_4": dict(sorted(Counter(r["det_mod_4"] for r in rows if r["even"]).items())),
               "even_odd_det_classes_by_type": dict(Counter(r["arf"]["type"] for r in rows if r["even"] and r["determinant"] % 2)),
               "decomposability_criteria_agree": True}
    return {"source": path, "source_sha256": hashlib.sha256(raw).hexdigest(),
            "universe": "the 169 isometry classes of positive definite integral rank-four lattices with det <= 24 "
                        "(accepted census, exact isometry dedup)", "lattices": rows, "summary": summary}


def task_a_eisenstein(path, max_nodes=5_000_000, progress=None):
    raw = _source(path).read_bytes()
    doc = json.loads(raw)
    reps, members = doc["representatives"], doc["members"]
    rows = []
    for n, r in enumerate(reps):
        m = members[r["member_index"]]
        G = m["L_reduced"]
        det = determinant(G)
        require(det == r["primitive_determinant"], "representative determinant")
        group = isometries(G, G, max_nodes)
        dec = decomposability(G, group, max_nodes)
        row = {"D_prime": r["D_prime"], "determinant": det, "det_mod_4": det % 4, "det_mod_8": det % 8,
               "even": is_even(G), "aut_order": len(group), "decomposable": dec["decomposable"],
               "eichler_components": dec["eichler"]["components"], "level": r["level"], "G": G, "a": m["a"], "b": m["b"],
               "member_index": r["member_index"], "P": "fails"}
        if is_even(G):
            row["arf"] = arf_type(G)
            row["f2"] = f2_data(G, m["a"], m["b"])
        rows.append(row)
        if progress and (n + 1) % 100 == 0:
            progress(f"  eisenstein {n + 1}/{len(reps)}")
    indec = [r for r in rows if not r["decomposable"]]
    smallest = {}
    for lvl in (1, 2):
        cands = [r for r in indec if r["level"] == lvl and not r["even"] and r["det_mod_4"] != 0]
        smallest[str(lvl)] = min(cands, key=lambda r: r["determinant"]) if cands else None
    return {"source": path, "source_sha256": hashlib.sha256(raw).hexdigest(),
            "universe": doc["universe"], "representatives": len(rows), "indecomposable": len(indec),
            "indecomposable_even_det_mod_4": dict(sorted(Counter(r["det_mod_4"] for r in indec if r["even"]).items())),
            "indecomposable_odd_lattice_det_mod_4": dict(sorted(Counter(r["det_mod_4"] for r in indec if not r["even"]).items())),
            "indecomposable_odd_lattice_det_not_div_4": [{k: r[k] for k in ("D_prime", "determinant", "aut_order", "level", "G", "a", "b")}
                                                         for r in sorted(indec, key=lambda r: r["determinant"])
                                                         if not r["even"] and r["det_mod_4"] != 0],
            "smallest_indecomposable_odd_lattice_det_not_div_4_by_level": smallest,
            "even_representatives_totally_singular": dict(Counter(str(r["f2"].get("span_a_c_totally_singular")) for r in rows if r["even"])),
            "lattices": rows}


def task_a_level4(pairs_path="results/rank4-level4/task-a.json",
                  certificate_paths=("examples/rank4-level4/det1729-level4.json",
                                     "examples/rank4-level4/other-level4-certificates.json"),
                  max_nodes=5_000_000, cutoff_vector_limit=5_000_000, progress=None):
    raw = _source(pairs_path).read_bytes()
    doc = json.loads(raw)
    members = []
    for pair in doc["pairs"]:
        for which, G in (("G1", pair["G1"]), ("G2", pair["G2"])):
            group = isometries(G, G, max_nodes)
            dec = decomposability(G, group, max_nodes)
            par = decide_parity(G, group, max_nodes, LADDER, cutoff_vector_limit)
            members.append({"pair": pair["name"], "which": which, "G": G, "determinant": determinant(G),
                            "det_mod_8": determinant(G) % 8, "even": is_even(G), "aut_order": len(group),
                            "decomposable": dec["decomposable"], "P": par["P"], "verdicts": par["verdicts"],
                            "largest_first_nonzero_norm4": par["largest_first_nonzero_norm4"],
                            "minimal_shell_decides": par["minimal_shell_decides"],
                            "totally_singular_pairs": sum(1 for d in par["decisions"] if d.get("f2", {}).get("span_a_c_totally_singular")),
                            "arf": par["arf"], "parity": par})
            if progress:
                progress(f"  pair member {pair['name']} {which}: P {par['P']} ({par['seconds']:.1f}s)")
    level4 = []
    for path in certificate_paths:
        craw = _source(path).read_bytes()
        cdoc = json.loads(craw)
        for cert in cdoc["certificates"]:
            G, a, b = cert["G"], cert["a"], cert["b"]
            group = isometries(G, G, max_nodes)
            dec = decomposability(G, group, max_nodes)
            par = decide_parity(G, group, max_nodes, LADDER, cutoff_vector_limit)
            known = next(d for d in par["decisions"] if d["a"] == a and d["b"] == b)
            if known["verdict"] == "unresolved":
                known["verdict"] = "zero_certified_externally"
                known["external_certificate"] = {"path": path, "sha256": hashlib.sha256(craw).hexdigest(), "level": cert["level"]}
                known["level"] = cert["level"]
                par["verdicts"] = dict(Counter(d["verdict"] for d in par["decisions"]))
                par["even_zeros"] = sum(1 for d in par["decisions"] if d["verdict"].startswith("zero"))
                par["zero_levels"] = dict(Counter(d.get("level") for d in par["decisions"] if d["verdict"].startswith("zero")))
                par["P"] = "fails" if par["verdicts"].get("unresolved", 0) == 0 or par["even_zeros"] else par["P"]
            level4.append({"source": path, "G": G, "a": a, "b": b, "determinant": cert["determinant"],
                           "det_mod_8": cert["determinant"] % 8, "even": is_even(G), "aut_order": len(group),
                           "decomposable": dec["decomposable"], "P": par["P"], "verdicts": par["verdicts"],
                           "zero_levels": par["zero_levels"], "arf": par["arf"],
                           "f2_of_level4_zero": f2_data(G, a, b) if is_even(G) else None, "parity": par})
            if progress:
                progress(f"  level-4 lattice det {cert['determinant']}: P {par['P']} verdicts {par['verdicts']} ({par['seconds']:.1f}s)")
    return {"pairs_source": pairs_path, "pairs_sha256": hashlib.sha256(raw).hexdigest(),
            "universe": "the eighteen members of the nine known isospectral pairs of RANK4-LEVEL4-SEARCH.md Task A "
                        "and the eight certified level-4 lattices", "pair_members": members, "level4_lattices": level4}


def task_a(max_nodes=5_000_000, cutoff_vector_limit=5_000_000, progress=None):
    start = time.perf_counter()
    census = task_a_census(max_nodes=max_nodes)
    if progress:
        progress("census done")
    eis_b = task_a_eisenstein("results/rank4-eisenstein/task-b.json", max_nodes, progress)
    eis_d = task_a_eisenstein("results/rank4-eisenstein/task-d.json", max_nodes, progress)
    lvl4 = task_a_level4(max_nodes=max_nodes, cutoff_vector_limit=cutoff_vector_limit, progress=progress)
    return {"schema": SCHEMA + "-task-a", "census": census, "eisenstein_task_b": eis_b, "eisenstein_task_d": eis_d,
            "level4": lvl4, "seconds": time.perf_counter() - start}


def render_task_a(doc):
    out = ["# Task A evidence table", "",
           "## 1. Census (169 classes, det <= 24): the 27 indecomposable classes", "",
           "| det | det mod 8 | parity | type | Aut | sqrt(-I) in Aut | P | even zeros | witness orders | first-nonzero norm4 max |",
           "|---:|---:|---|---|---:|---|---|---:|---|---:|"]
    for r in sorted((r for r in doc["census"]["lattices"] if not r["decomposable"]), key=lambda r: (r["determinant"], r["even"])):
        out.append(f"| {r['determinant']} | {r['det_mod_8']} | {'even' if r['even'] else 'odd'} | {(r.get('arf') or {}).get('type', '-')} | "
                   f"{r['aut_order']} | {r['has_square_root_of_minus_identity']} | {r['P']} | {r['even_zeros']} | "
                   f"{r['zero_witness_orders']} | {max(map(int, r['first_nonzero_norm4_distribution']), default='-')} |")
    s = doc["census"]["summary"]
    out += ["", f"Decomposable classes: {s['classes'] - s['indecomposable']}, all fail P: {s['every_decomposable_class_fails_P']}.",
            f"Holding: {s['indecomposable_holding']}.  Failing: {s['indecomposable_failing']}.",
            f"Even classes by det mod 4: {s['even_classes_det_mod_4']}.  Even odd-det classes by type: {s['even_odd_det_classes_by_type']}.", ""]
    for key, title in (("eisenstein_task_b", "Eisenstein Task B representatives"), ("eisenstein_task_d", "Eisenstein Task D representatives")):
        e = doc[key]
        out += [f"## {title} ({e['representatives']} classes, all fail P)", "",
                f"Indecomposable: {e['indecomposable']}; even ones by det mod 4: {e['indecomposable_even_det_mod_4']}; "
                f"odd-lattice ones by det mod 4: {e['indecomposable_odd_lattice_det_mod_4']}.", "",
                "Indecomposable odd lattices with 4 not dividing det:", "",
                "| D' | det | Aut | level | G | a | b |", "|---:|---:|---:|---:|---|---|---|"]
        for r in e["indecomposable_odd_lattice_det_not_div_4"]:
            out.append(f"| {r['D_prime']} | {r['determinant']} | {r['aut_order']} | {r['level']} | {r['G']} | {r['a']} | {r['b']} |")
        sm = e["smallest_indecomposable_odd_lattice_det_not_div_4_by_level"]
        out += ["", "Smallest by level: " + "; ".join(f"level {k}: det {v['determinant'] if v else None}" for k, v in sm.items()), ""]
    out += ["## Isospectral pair members and level-4 lattices", "",
            "| lattice | det | det mod 8 | parity | type | Aut | decomposable | P | verdicts | largest first-nonzero norm4 | ts pairs |",
            "|---|---:|---:|---|---|---:|---|---|---|---:|---:|"]
    for r in doc["level4"]["pair_members"]:
        out.append(f"| {r['pair']} {r['which']} | {r['determinant']} | {r['det_mod_8']} | {'even' if r['even'] else 'odd'} | "
                   f"{(r['arf'] or {}).get('type', '-')} | {r['aut_order']} | {r['decomposable']} | {r['P']} | {r['verdicts']} | "
                   f"{r['largest_first_nonzero_norm4']} | {r['totally_singular_pairs']} |")
    for r in doc["level4"]["level4_lattices"]:
        out.append(f"| level-4 det {r['determinant']} | {r['determinant']} | {r['det_mod_8']} | {'even' if r['even'] else 'odd'} | "
                   f"{(r['arf'] or {}).get('type', '-')} | {r['aut_order']} | {r['decomposable']} | {r['P']} | {r['verdicts']} | "
                   f"{r['parity']['largest_first_nonzero_norm4']} | - |")
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------------
# Task D: structure, assembled from the data
# ----------------------------------------------------------------------------

def complex_structure_plane(G, group):
    """For g with g^2 = -I: the plane (g+I)L/2L, which is totally singular."""
    for g in group:
        if multiply(g, g) == MINUS_I4:
            cols = transpose([[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(g, I4)])
            classes = sorted({tuple(sum(k * c[i] for k, c in zip(coeffs, cols)) % 2 for i in range(4))
                              for coeffs in product((0, 1), repeat=4)} - {(0, 0, 0, 0)})
            return {"g": g, "plane_nonzero_classes": [list(c) for c in classes],
                    "totally_singular": all(quadratic_value(G, list(c)) == 0 for c in classes)}
    return None


def task_d(task_a_doc, task_b_doc, task_c_doc=None):
    zeros = []

    def add(source, det, decomposable, level, f2, extra=None):
        zeros.append({"source": source, "determinant": det, "det_mod_8": det % 8, "decomposable": decomposable,
                      "level": level, "type": f2["type"], "radical_dimension": f2["radical_dimension"],
                      "radical_totally_singular": f2["radical_totally_singular"], "a_in_radical": f2["a_in_radical"],
                      "b_in_image": f2["b_in_image_mod_2"], "Q_a": f2["Q_a"], "Q_c": f2.get("Q_c_values"),
                      "a_equals_c": f2.get("a_equals_some_c"), "span_totally_singular": f2.get("span_a_c_totally_singular"),
                      "a_in_radical_plus_c": f2.get("a_in_radical_plus_c"), **(extra or {})})
    for r in task_a_doc["census"]["lattices"]:
        if r["even"]:
            for z in r["zeros"]:
                add("census", r["determinant"], r["decomposable"], 1, z["f2"], {"witness_order": z["witness_order"]})
    for key in ("eisenstein_task_b", "eisenstein_task_d"):
        for r in task_a_doc[key]["lattices"]:
            if r["even"]:
                add(key, r["determinant"], r["decomposable"], r["level"], r["f2"])
    for r in task_a_doc["level4"]["level4_lattices"]:
        if r["even"]:
            for d in r["parity"]["decisions"]:
                if d["verdict"].startswith("zero"):
                    add("level4", r["determinant"], r["decomposable"], d.get("level"), d["f2"])
    for doc, name in ((task_b_doc, "task_b"), (task_c_doc, "task_c")):
        if doc is None:
            continue
        for r in doc["lattices"]:
            if "parity" in r:
                for d in r["parity"]["decisions"]:
                    if d["verdict"].startswith("zero"):
                        add(name, r["determinant"], r["decomposable"], d.get("level"), d["f2"],
                            {"witness_order": d.get("witness_order")})
    indec_zeros = [z for z in zeros if not z["decomposable"]]
    d1 = {"known_even_zeros_of_even_lattices": len(zeros), "indecomposable": len(indec_zeros),
          "indecomposable_all_span_totally_singular": all(z["span_totally_singular"] for z in indec_zeros),
          "indecomposable_violations": [z for z in indec_zeros if not z["span_totally_singular"]],
          "decomposable_violations": sum(1 for z in zeros if z["decomposable"] and not z["span_totally_singular"]),
          "by_source_and_det_mod_8": dict(Counter(f"{z['source']} det%8={z['det_mod_8']} level {z['level']}" for z in zeros)),
          "degenerate_indecomposable_positions": dict(Counter(
              f"rad dim {z['radical_dimension']}, rad ts {z['radical_totally_singular']}, a in rad {z['a_in_radical']}, "
              f"b in im {z['b_in_image']}, Q(a) {z['Q_a']}, Q(c) {z['Q_c']}, ts {z['span_totally_singular']}"
              for z in indec_zeros if z["type"] == "degenerate"))}
    # D2: plus-type classes with odd det
    plus = []
    sources = [("census", r) for r in task_a_doc["census"]["lattices"] if r["even"] and r["determinant"] % 2 and not r["decomposable"]]
    sources += [("task_b", r) for r in task_b_doc["lattices"] if not r["decomposable"]]
    sources += [("pair_member", r) for r in task_a_doc["level4"]["pair_members"] if r["even"] and not r["decomposable"]]
    for name, r in sources:
        if name == "census":
            ts, zs, has_sqrt, orders = r["totally_singular_pairs"], r["even_zeros"], r["has_square_root_of_minus_identity"], r["zero_witness_orders"]
        else:
            par = r["parity"]
            ts = sum(1 for d in par["decisions"] if d.get("f2", {}).get("span_a_c_totally_singular"))
            zs = par["even_zeros"]
            has_sqrt = r.get("has_square_root_of_minus_identity")
            if has_sqrt is None and name == "pair_member":
                has_sqrt = None
            orders = dict(Counter(d.get("witness_order") for d in par["decisions"] if d["verdict"] == "zero_by_symmetry"))
        plus.append({"source": name, "determinant": r["determinant"], "det_mod_8": r["det_mod_8"],
                     "type": (r["arf"] or {}).get("type"), "aut_order": r["aut_order"],
                     "has_square_root_of_minus_identity": has_sqrt, "totally_singular_ordered_pairs": ts,
                     "even_zeros": zs, "zero_witness_orders": orders, "P": r["P"] if name != "task_b" else r["parity"]["P"]})
    d2 = {"classes": plus,
          "plus_type_classes": sum(1 for p in plus if p["type"] == "plus"),
          "plus_type_with_zeros": sum(1 for p in plus if p["type"] == "plus" and p["even_zeros"]),
          "plus_type_zero_classes_all_have_sqrt_minus_I": all(p["has_square_root_of_minus_identity"] for p in plus if p["type"] == "plus" and p["even_zeros"]),
          "plus_type_with_sqrt_minus_I_all_fail": all(p["even_zeros"] == 6 for p in plus if p["type"] == "plus" and p["has_square_root_of_minus_identity"]),
          "minus_type_classes": sum(1 for p in plus if p["type"] == "minus"),
          "minus_type_with_zeros": sum(1 for p in plus if p["type"] == "minus" and p["even_zeros"]),
          "minus_type_with_sqrt_minus_I": sum(1 for p in plus if p["type"] == "minus" and p["has_square_root_of_minus_identity"]),
          "totally_singular_pairs_per_plus_class": dict(Counter(p["totally_singular_ordered_pairs"] for p in plus if p["type"] == "plus")),
          "involution_witnesses_in_indecomposable_odd_det": sum(v for p in plus for k, v in p["zero_witness_orders"].items() if k == 2)}
    d2["classes_with_sqrt_minus_I"] = [{"determinant": r["determinant"], "aut_order": r["aut_order"], "decomposable": r["decomposable"],
                                        "P": r["parity"]["P"] if "parity" in r else "fails (decomposable)"}
                                       for r in task_b_doc["lattices"] if r["has_square_root_of_minus_identity"]]
    d2["every_class_with_sqrt_minus_I_fails_P"] = all(c["P"].startswith("fails") for c in d2["classes_with_sqrt_minus_I"])
    d3 = []
    for r in task_b_doc["lattices"]:
        if "parity" in r and r["parity"]["P"] == "fails":
            plane = complex_structure_plane(r["G"], r["automorphisms"])
            zs = [d for d in r["parity"]["decisions"] if d["verdict"].startswith("zero")]
            spans = sorted({tuple(d["a"]) for d in zs} | {tuple(c) for d in zs for c in d["f2"]["c_solutions"]})
            d3.append({"determinant": r["determinant"], "G": r["G"], "aut_order": r["aut_order"], "complex_structure": plane,
                       "zero_classes_a_and_c": [list(s) for s in spans],
                       "zeros_lie_on_the_plane": plane is not None and sorted(map(tuple, plane["plane_nonzero_classes"])) == spans,
                       "zero_count": len(zs), "witness_orders": dict(Counter(d.get("witness_order") for d in zs)),
                       "every_witness_squares_to_minus_identity": all(multiply(d["witness"], d["witness"]) == MINUS_I4 for d in zs if "witness" in d),
                       "zeros": [{"a": d["a"], "b": d["b"], "c": d["f2"]["c_solutions"], "level": d.get("level")} for d in zs]})
    return {"schema": SCHEMA + "-task-d", "d1": d1, "d2": d2, "d3": d3, "zeros": zeros}


# ----------------------------------------------------------------------------
# Task D3: the E8 witness is a Z[i]-structure (checked example)
# ----------------------------------------------------------------------------

def integer_basis(generators):
    """Column basis of the Z-span of integer vectors (generic Hermite reduction)."""
    cols = [list(g) for g in generators]
    n = len(cols[0])
    basis = []
    for i in range(n):
        pivots = [c for c in cols if c[i] != 0]
        rest = [c for c in cols if c[i] == 0]
        while len(pivots) > 1:
            pivots.sort(key=lambda c: abs(c[i]))
            p = pivots[0]
            new = [p]
            for c in pivots[1:]:
                q = c[i] // p[i]
                c = [x - q * y for x, y in zip(c, p)]
                (new if c[i] != 0 else rest).append(c)
            pivots = new
        if pivots:
            p = pivots[0]
            basis.append([-x for x in p] if p[i] < 0 else p)
        cols = rest
    return transpose(basis)


def e8_check():
    """E8 = D8^+ carries g = four coordinate rotations, g^2 = -I; its fixed four-space is totally singular
    and every pair of distinct nonzero classes in it is a symmetry zero with witness g."""
    n = 8
    gens = [[2 * (int(k == i) - int(k == i + 1)) for k in range(n)] for i in range(n - 1)]
    gens.append([2 * (int(k == n - 2) + int(k == n - 1)) for k in range(n)])
    gens.append([1] * n)                                    # 2 * (1/2, ..., 1/2)
    B2 = integer_basis(gens)                                # columns: 2 * (basis of E8)
    G4 = multiply(transpose(B2), B2)
    require(all(x % 4 == 0 for row in G4 for x in row), "E8 Gram")
    G = [[x // 4 for x in row] for row in G4]
    require(determinant(G) == 1 and all(G[i][i] % 2 == 0 for i in range(n)), "E8 is even unimodular")
    J = [[0] * n for _ in range(n)]
    for k in range(0, n, 2):
        J[k][k + 1], J[k + 1][k] = -1, 1                      # (x, y) -> (-y, x) on each coordinate pair
    g = multiply(inverse(B2), multiply(J, B2))
    require(all(Fraction(x).denominator == 1 for row in g for x in row), "rotation preserves E8")
    g = [[int(Fraction(x)) for x in row] for row in g]
    I = identity(n)
    require(multiply(multiply(transpose(g), G), g) == G and multiply(g, g) == [[-x for x in row] for row in I], "g^2 = -I")
    cols = transpose([[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(g, I)])
    classes = sorted({tuple(sum(k * c[i] for k, c in zip(coeffs, cols)) % 2 for i in range(n))
                      for coeffs in product((0, 1), repeat=n)})
    nonzero = [list(c) for c in classes if any(c)]
    require(len(classes) == 16, "fixed space has dimension 4")
    singular = all(dot(c, matvec(G, c)) % 4 == 0 for c in nonzero)
    witnessed, orthogonal = 0, 0
    for a in nonzero:
        for c in nonzero:
            if a == c:
                continue
            b = [x % 2 for x in matvec(G, c)]
            require(dot(a, b) % 2 == 0, "fixed space not totally singular")
            orthogonal += 1
            # epsilon(g) = B(alpha, c) with a = (g+I) alpha: in rank 8 this is 1 for 8 of the 14 partners of c
            witnessed += symmetry_data(G, a, b, g)["epsilon"] == 1
    checked = orthogonal
    seen = set()
    for a in nonzero:
        for c in nonzero:
            if a != c:
                seen.add(frozenset((tuple(a), tuple(c), tuple((x + y) % 2 for x, y in zip(a, c)))))
    planes = len(seen)
    return {"E8_gram": G, "g": g, "g_squared_is_minus_identity": True, "fixed_space_dimension": 4,
            "fixed_space_totally_singular": singular, "planes_in_fixed_space": planes,
            "ordered_pairs_in_fixed_space": checked, "ordered_pairs_with_epsilon_one_for_g": witnessed,
            "aut_E8_order": 696729600,
            "note": "every ordered pair of distinct nonzero classes of the fixed four-space is an E8 zero by the plane criterion "
                    "of E7-E8-MILESTONE.md; g itself is the epsilon-one witness for 8 of the 14 partners of each class, the "
                    "others need conjugates of g (all 1575 planes form one Aut(E8)-orbit).  In rank four the fixed space is a "
                    "plane and g witnesses all six ordered pairs."}


# ----------------------------------------------------------------------------
# command line
# ----------------------------------------------------------------------------

def fresh(path):
    path = Path(path)
    require(not path.exists(), f"choose a fresh output path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def save(path, doc):
    fresh(path).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("A", "B", "C", "D", "E", "checks", "certify", "estimate", "e8"))
    parser.add_argument("--dmax", type=int, default=49)
    parser.add_argument("--dmin", type=int, default=1)
    parser.add_argument("--max-nodes", type=int, default=5_000_000)
    parser.add_argument("--cutoff-limit", type=int, default=5_000_000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--task-a", type=Path, default=Path("results/rank4-parity/task-a.json"))
    parser.add_argument("--task-b", type=Path)
    parser.add_argument("--task-c", type=Path)
    parser.add_argument("--examples-dir", type=Path, default=Path("examples/rank4-parity"))
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    import sys
    progress = (lambda msg: print(msg, file=sys.stderr, flush=True)) if args.progress else None
    if args.task == "checks":
        doc = generator_checks(max_nodes=args.max_nodes)
    elif args.task == "A":
        doc = task_a(args.max_nodes, args.cutoff_limit, progress)
        if args.markdown:
            fresh(args.markdown).write_text(render_task_a(doc), encoding="utf-8")
    elif args.task == "B":
        doc = task_b(args.dmax, args.max_nodes, args.cutoff_limit, progress)
    elif args.task == "C":
        doc = task_c(args.dmax, args.max_nodes, args.cutoff_limit, progress)
    elif args.task == "E":
        doc = task_e(args.dmin, args.dmax, args.max_nodes, args.cutoff_limit, progress)
    elif args.task == "estimate":
        doc = {"schema": SCHEMA + "-task-e-estimate", "rows": task_e_estimate()}
    elif args.task == "e8":
        doc = {"schema": SCHEMA + "-e8-check", **e8_check()}
    elif args.task == "D":
        a = json.loads(args.task_a.read_text(encoding="utf-8"))
        b = json.loads(args.task_b.read_text(encoding="utf-8"))
        c = json.loads(args.task_c.read_text(encoding="utf-8")) if args.task_c else None
        doc = task_d(a, b, c)
    elif args.task == "certify":
        b = json.loads(args.task_b.read_text(encoding="utf-8"))
        certs = []
        for r in b["lattices"]:
            if "parity" in r and r["parity"]["P"] == "fails":
                for d in r["parity"]["decisions"]:
                    if d["verdict"].startswith("zero"):
                        certs.append(certify_counterexample(r["G"], d["a"], d["b"], r["automorphisms"], args.max_nodes))
                        if progress:
                            progress(f"  certified det {r['determinant']} a {d['a']} b {d['b']} level {certs[-1]['level']}")
        doc = {"schema": SCHEMA + "-c1-counterexamples", "source_task_b": str(args.task_b), "certificates": certs}
    if args.output:
        save(args.output, doc)
    summary = doc.get("summary", {k: v for k, v in doc.items() if k not in ("lattices", "zeros", "certificates", "classes_doc", "census", "eisenstein_task_b", "eisenstein_task_d", "level4", "candidate_bases", "representatives", "generation")})
    if args.task == "A":
        summary = {"census": doc["census"]["summary"], "eisenstein_b": doc["eisenstein_task_b"]["smallest_indecomposable_odd_lattice_det_not_div_4_by_level"],
                   "eisenstein_d": doc["eisenstein_task_d"]["smallest_indecomposable_odd_lattice_det_not_div_4_by_level"],
                   "pair_members_P": dict(Counter(m["P"] for m in doc["level4"]["pair_members"])),
                   "level4_P": dict(Counter(m["P"] for m in doc["level4"]["level4_lattices"])), "seconds": doc["seconds"]}
    if args.task == "certify":
        summary = {"certificates": [(c["determinant"], c["a"], c["b"], c["level"]) for c in doc["certificates"]]}
    print(json.dumps(summary, default=str))


if __name__ == "__main__":
    main()
