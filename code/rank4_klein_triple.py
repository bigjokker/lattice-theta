"""Exact check of the Klein-configuration reading of the rank-four counterexample.

For an even characteristic (a,b) on a lattice L with Gram matrix G, let
L0 = {x : x.b even}, t in L \ L0, xi = a/2, L1 = L0 + Z xi, L2 = L0 + Z(xi + t)
and L' = L + Z xi.  Then Theta_G[a,b] = i^(a.b) (theta_L1 - theta_L2), so the
characteristic vanishes identically exactly when L1 and L2 are isospectral.

This script verifies, with exact rational arithmetic and no third-party
dependencies, that in the determinant-23,716 counterexample the three
intermediate lattices L, L1, L2 between L0 and L' are pairwise isometric and
that the ambient order-three rotation R preserves L0 and permutes them
cyclically: L1 -> L -> L2 -> L1.  It prints their Gram matrices.

Run:  python rank4_klein_triple.py
Exit status 0 means every assertion held.
"""
from fractions import Fraction as Fr

G = [[12, 12, -4, -1], [12, 16, -6, 2], [-4, -6, 28, -14], [-1, 2, -14, 28]]
M = [[8, -4, 1, 2], [-4, 8, -3, 1], [1, -3, 14, -7], [2, 1, -7, 14]]
V = [[-1, 0, 0, 0], [1, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]
R = [[0, -1, 0, 0], [1, -1, 0, 0], [0, 0, 0, -1], [0, 0, 1, -1]]
A = [0, 1, 0, 0]   # a
B = [1, 0, 0, 0]   # b, so L0 = {x : x0 even} and t = e1


def mul(X, Y):
    return [[sum(Fr(X[i][k]) * Y[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def tr(X):
    return [list(r) for r in zip(*X)]


def inv(X):
    n = 4
    W = [[Fr(v) for v in r] + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(X)]
    for i in range(n):
        p = next(r for r in range(i, n) if W[r][i] != 0)
        W[i], W[p] = W[p], W[i]
        W[i] = [v / W[i][i] for v in W[i]]
        for r in range(n):
            if r != i and W[r][i] != 0:
                W[r] = [W[r][c] - W[r][i] * W[i][c] for c in range(2 * n)]
    return [r[n:] for r in W]


def det(X):
    W = [[Fr(v) for v in r] for r in X]
    d = Fr(1)
    for i in range(4):
        p = next((r for r in range(i, 4) if W[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            W[i], W[p] = W[p], W[i]
            d = -d
        d *= W[i][i]
        for r in range(i + 1, 4):
            f = W[r][i] / W[i][i]
            W[r] = [W[r][c] - f * W[i][c] for c in range(4)]
    return d


def integral(X):
    return all(v.denominator == 1 for r in X for v in r)


def columns(*cols):
    return [list(r) for r in zip(*cols)]


def main():
    assert sum(x * y for x, y in zip(A, B)) % 2 == 0, "characteristic must be even"
    assert mul(mul(tr(V), M), V) == [[2 * v for v in r] for r in G]
    assert mul(mul(tr(R), M), R) == [[Fr(v) for v in r] for r in M]

    # Bases of the doubled lattices 2X in the coordinates y = 2x + a.
    # 2L0 = diag(4,2,2,2) Z^4, t = e1, a = (0,1,0,0).
    doubled = {
        "L0": columns((4, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)),
        "L": columns((2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)),
        "L1": columns((4, 0, 0, 0), (0, 1, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)),
        "L2": columns((4, 0, 0, 0), (2, 1, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)),
    }

    def gram(name):
        Bm = doubled[name]
        return [[v / 4 for v in r] for r in mul(mul(tr(Bm), G), Bm)]

    for name in ("L", "L1", "L2", "L0"):
        g = gram(name)
        print(name, "Gram", [[int(v) for v in r] for r in g], "det", det(g))
        assert integral(g)
    assert gram("L") == [[Fr(v) for v in r] for r in G]
    assert det(gram("L")) == det(gram("L1")) == det(gram("L2")) == 23716
    assert det(gram("L0")) == 4 * 23716

    def r_maps(src, dst):
        # R acts on ambient coordinates w = V y; check R(V B_src Z^4) = V B_dst Z^4.
        U = mul(inv(mul(V, doubled[dst])), mul(R, mul(V, doubled[src])))
        return integral(U) and abs(det(U)) == 1

    assert r_maps("L0", "L0"), "R must preserve L0"
    assert r_maps("L1", "L"), "R must map L1 onto L"
    assert r_maps("L", "L2"), "R must map L onto L2"
    assert r_maps("L2", "L1"), "R must map L2 onto L1"
    print("R preserves L0 and permutes L1 -> L -> L2 -> L1.")

    # Explicit unimodular U with U^t G U = Gram(L1): coordinates of R(basis of L1)
    # in the basis of L, read back in the original coordinates.
    U = mul(inv(mul(V, doubled["L"])), mul(R, mul(V, doubled["L1"])))
    assert integral(U) and abs(det(U)) == 1
    assert mul(mul(tr(U), G), U) == gram("L1")
    print("Unimodular U with U^t G U = Gram(L1):", [[int(v) for v in r] for r in U])
    print("All assertions held: L, L1 and L2 are pairwise isometric.")


if __name__ == "__main__":
    main()
