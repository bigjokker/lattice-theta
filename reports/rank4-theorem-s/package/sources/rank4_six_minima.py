"""Six sub-coset minima in the common lattice K, and exact symmetry certificates (RANK4-SIX-MINIMA-COMMON-LATTICE.md).

Setting: L even, positive definite, rank four, odd determinant, Gram G; P a B-isotropic plane of L/2L
given by its three nonzero class indices (BITS labels, as in rank4_pairing_search).

  K = {y in L : y mod 2L in P}, index 4, containing 2L; basis columns B (integer HNF of 2L + a + c).
  GK = B^t G B (all entries even), H = GK / 2, det H = det G.
  The 16 classes of K/2K: the four classes inside W = 2L/2K are fixed; the other twelve form six pairs,
  one per (shift a in P, half eps = B(a, x) mod 2), and the two classes of a pair carry opposite signs.
  tau = the involution of K/2K swapping each pair (linear, tau - I has image and kernel W).

For each of the six (shift, half): the signed half-series equals theta(v + 2K) - theta(w + 2K) at every
coefficient, so its minimal shell is balanced iff the two ordinary cosets have equal minimum AND equal
minimal multiplicity.  six_minima() returns the six records with complete minimizing sets; every minimum
is reached (exact radius: an explicit small vector of each class bounds that class's minimum).

Certificates.
  marked zero:  an integral T with T^t GK T = GK and T mod 2 = tau proves all six coset identities, hence
                the three signed theta series of the plane vanish identically (all orders).
  splitting:    an integral h != +-I with h^2 = I, h = I mod 2, h^t G h = G; then (I+h)/2, (I-h)/2 are
                integral orthogonal projections, so L is decomposable.
check_marked_certificate() and check_splitting() verify a stored certificate from (G, plane, B, T) alone.

Usage
  python -B rank4_six_minima.py profile --G "[[6,0,1,2],[0,6,-2,1],[1,-2,16,0],[2,1,0,18]]" [--plane 3 12 15]
  python -B rank4_six_minima.py certify --G "[[...]]" --plane 5 10 15
"""
import argparse
from itertools import product
import json

from exact_theta import InvalidCertificate, dot, inverse, matvec, multiply, require, transpose
import rank4_pairing_search as S

BITS = S.BITS


def label(v):
    return sum((int(x) % 2) << i for i, x in enumerate(v))


def _int_matrix(M, name):
    require(isinstance(M, list) and len(M) == 4 and all(isinstance(r, list) and len(r) == 4 for r in M),
            f"{name} must be a 4x4 list of lists")
    require(all(type(x) is int for r in M for x in r), f"{name} must have exact integer entries")
    return M


def validate_gram(G, even_odd=True):
    """G must be an integral symmetric positive definite 4x4 matrix; with even_odd also even of odd determinant."""
    _int_matrix(G, "G")
    require(all(G[i][j] == G[j][i] for i in range(4) for j in range(4)), "G must be symmetric")
    require(all(S.determinant([r[:k] for r in G[:k]]) > 0 for k in range(1, 5)), "G must be positive definite")
    if even_odd:
        require(all(G[i][i] % 2 == 0 for i in range(4)) and S.determinant(G) % 2 == 1, "G must be even of odd determinant")
    return G


def _check_plane(G, plane):
    plane = tuple(sorted(int(k) for k in plane))
    require(len(plane) == 3 and len(set(plane)) == 3 and all(0 < k < 16 for k in plane), "plane: three nonzero classes")
    require(plane[0] ^ plane[1] == plane[2], "plane: not a 2-dimensional subspace")
    require(all(sum(BITS[u][i] * G[i][j] * BITS[v][j] for i in range(4) for j in range(4)) % 2 == 0
                for u in plane for v in plane), "plane: not B-isotropic")
    return plane


def common_lattice(G, plane):
    """(plane, B, GK, H, W, pairs, tau).  pairs: {(shift, half): [(class, sign), (class, sign)]} with signs +1, -1."""
    validate_gram(G)
    plane = _check_plane(G, plane)
    a, c = [list(BITS[k]) for k in plane[:2]]
    gens = [[2 * int(i == j) for i in range(4)] for j in range(4)] + [a, c]
    B = transpose(S._int_hnf(gens))
    require(abs(S.determinant(B)) == 4, "K must have index 4")
    GK = multiply(multiply(transpose(B), G), B)
    require(all(x % 2 == 0 for r in GK for x in r), "B^t G B must be even")
    H = [[x // 2 for x in r] for r in GK]
    require(S.determinant(H) == S.determinant(G), "det H must equal det G")
    W, pairs = [], {}
    for k in range(16):
        y = matvec(B, list(BITS[k]))
        ak = label(y)
        if ak == 0:
            W.append(k)
            continue
        require(ak in plane, "a K-class outside the plane")
        av, cv = list(BITS[ak]), list(BITS[next(j for j in plane if j != ak)])
        x = [(p - q) // 2 for p, q in zip(y, av)]
        half = dot(av, matvec(G, x)) % 2
        sign = 1 - 2 * (dot(cv, matvec(G, x)) % 2)
        pairs.setdefault((ak, half), []).append((k, sign))
    require(len(W) == 4 and len(pairs) == 6, "K/2K must have 4 fixed classes and 6 pairs")
    tau = list(range(16))
    for key, pr in pairs.items():
        require(len(pr) == 2 and pr[0][1] == -pr[1][1], "each pair must carry opposite signs")
        pr.sort(key=lambda t: -t[1])                 # plus class first
        p, q = pr[0][0], pr[1][0]
        tau[p], tau[q] = q, p
    require(all(tau[p ^ q] == tau[p] ^ tau[q] for p in range(16) for q in range(16)), "tau must be linear")
    require(all(tau[tau[p]] == p for p in range(16)) and {p for p in range(16) if tau[p] == p} == set(W),
            "tau must be an involution fixing exactly W")
    return plane, B, GK, H, W, pairs, tau


def class_minimum(GK, k, max_nodes=50_000_000):
    """(minimum norm, all minimal vectors) of the coset BITS[k] + 2K, in K coordinates (norms are L-norms)."""
    v0 = list(BITS[k])
    R = None
    for i in range(-1, 4):                            # explicit small vectors bound the minimum
        for j in range(i, 4):
            for si, sj in product((-1, 1), repeat=2):
                z = v0[:]
                if i >= 0:
                    z[i] += 2 * si
                if j >= 0 and j != i:
                    z[j] += 2 * sj
                if any(z):
                    n = dot(z, matvec(GK, z))
                    R = n if R is None else min(R, n)
    found = []

    def inner(lo, hi, step, cc, C, y):
        for y0 in range(lo, hi + 1, step):
            z = [y0, y[1], y[2], y[3]]
            if any(z):
                found.append((dot(z, matvec(GK, z)), tuple(z)))
    S._enumerate(GK, R, v0, inner, max_nodes)
    require(found, "class minimum: no vector found within the bound")
    m = min(n for n, _ in found)
    return m, sorted(z for n, z in found if n == m)


def six_minima(G, plane):
    """The six (shift, half) records of the plane, each with both coset minima and complete minimizing sets."""
    plane, B, GK, H, W, pairs, tau = common_lattice(G, plane)
    out = []
    for (ak, half), ((kp, _), (km, _)) in sorted(pairs.items()):
        mp, vp = class_minimum(GK, kp)
        mm, vm = class_minimum(GK, km)
        m = min(mp, mm)
        plus = len(vp) if mp == m else 0
        minus = len(vm) if mm == m else 0
        out.append({"shift": ak, "half": half, "residue_mod_8": (S.norm(G, list(BITS[ak])) + 4 * half) % 8,
                    "min_norm4": m, "plus": plus, "minus": minus, "signed": plus - minus,
                    "balanced": mp == mm and len(vp) == len(vm), "reached": True,
                    "plus_class": kp, "minus_class": km, "plus_coset_min": mp, "minus_coset_min": mm,
                    "plus_minimizers_L": [matvec(B, list(z)) for z in vp] if mp == m else [],
                    "minus_minimizers_L": [matvec(B, list(z)) for z in vm] if mm == m else []})
    return out


def find_marked_certificate(G, plane, max_nodes=20_000_000):
    """An automorphism T of K lifting tau, as a certificate {G, plane, B, T, h}; None if Aut(K) has no lift."""
    from verify_rank4_feasibility import box_isometries
    plane, B, GK, H, W, pairs, tau = common_lattice(G, plane)
    for T in box_isometries(H, H, max_nodes):
        if S.mod2_perm(T) == tau:
            h = multiply(multiply(B, T), inverse(B))
            return {"G": G, "plane": list(plane), "B": B, "T": T, "h": [[int(x) for x in r] for r in h]}
    return None


def check_marked_certificate(cert):
    """Independent check: B spans K, T is an integral automorphism of B^t G B, and T mod 2 is tau.
    Returns True or raises InvalidCertificate.  Implies the three signed series of the plane vanish identically."""
    G, B, T = validate_gram(cert["G"]), _int_matrix(cert["B"], "B"), _int_matrix(cert["T"], "T")
    plane = _check_plane(G, cert["plane"])
    require(abs(S.determinant(B)) == 4, "certificate: [L:K] must be 4")
    require(all(label([B[i][j] for i in range(4)]) in (0,) + plane for j in range(4)), "certificate: B outside K")
    GK = multiply(multiply(transpose(B), G), B)
    require(abs(S.determinant(T)) == 1 and multiply(multiply(transpose(T), GK), T) == GK, "certificate: T not in Aut(K)")
    _, B2, _, _, _, _, tau = common_lattice(G, plane)
    # tau in the coordinates of the supplied basis B: conjugate by the change of basis B2 -> B (integral, det +-1)
    C = multiply(inverse(B2), B)
    require(all(x.denominator == 1 for r in C for x in r) and abs(S.determinant([[int(x) for x in r] for r in C])) == 1,
            "certificate: B is not a basis of K")
    Cint = [[int(x) for x in r] for r in C]                    # K coordinates: B-coords -> B2-coords
    perm_T = S.mod2_perm(multiply(multiply(Cint, T), [[int(x) for x in r] for r in inverse(Cint)]))
    require(perm_T == tau, "certificate: T mod 2 is not tau")
    h = multiply(multiply(B, T), inverse(B))
    require(all(x.denominator == 1 for r in h for x in r), "certificate: T does not extend to L")
    hi = [[int(x) for x in r] for r in h]
    require(multiply(multiply(transpose(hi), G), hi) == G, "certificate: extension is not an isometry of L")
    return True


def check_splitting(G, h):
    """Integral splitting witness: h != +-I, h^2 = I, h = I mod 2, h^t G h = G.  True or InvalidCertificate."""
    validate_gram(G, even_odd=False)
    _int_matrix(h, "h")
    I = [[int(i == j) for j in range(4)] for i in range(4)]
    require(h != I and h != [[-x for x in r] for r in I], "splitting: h must not be +-I")
    require(multiply(h, h) == I, "splitting: h^2 must be I")
    require(all((h[i][j] - I[i][j]) % 2 == 0 for i in range(4) for j in range(4)), "splitting: h must be I mod 2")
    require(multiply(multiply(transpose(h), G), h) == G, "splitting: h must be an isometry")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", choices=("profile", "certify"))
    ap.add_argument("--G", required=True, help="Gram matrix as JSON")
    ap.add_argument("--plane", type=int, nargs=3)
    a = ap.parse_args()
    G = json.loads(a.G)
    planes = [tuple(a.plane)] if a.plane else [tuple(p) for p in S.isotropic_planes(G)]
    for p in planes:
        rows = six_minima(G, p)
        rec = {"plane": list(p), "balanced": sum(r["balanced"] for r in rows),
               "six": [{k: r[k] for k in ("shift", "half", "residue_mod_8", "min_norm4", "plus", "minus", "balanced")}
                       for r in rows]}
        if a.task == "certify" and rec["balanced"] == 6:
            cert = find_marked_certificate(G, p)
            rec["certificate"] = cert
            rec["certificate_checked"] = cert is not None and check_marked_certificate(cert)
        print(json.dumps(rec))


if __name__ == "__main__":
    main()
