"""Exact rank-three splitting and obtuse-superbase shell certificates."""

from itertools import combinations, product
from math import gcd

from exact_theta import (EnumerationLimit, dot, identity, inverse, ldl, matrix,
                         integer, matvec, multiply, require, transpose, vector)
from lattice_census import sphere


def admit(G):
    matrix(G, "G", 3)
    ldl(G)


def egcd(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def complete_primitive(v):
    """Integer row Bezout operations send v to e1; invert to complete v."""
    vector(v, "primitive vector", 3)
    P, w = identity(3), list(v)
    for i in (1, 2):
        a, b = w[0], w[i]
        if not a and not b:
            continue
        g, s, t = egcd(a, b)
        row0, rowi = P[0][:], P[i][:]
        P[0] = [s*x + t*y for x, y in zip(row0, rowi)]
        P[i] = [-b//g*x + a//g*y for x, y in zip(row0, rowi)]
        w[0], w[i] = g, 0
    require(w == [1, 0, 0], "v must be primitive")
    U = inverse(P)
    require(all(x.denominator == 1 for row in U for x in row), "completion failed")
    return [[int(x) for x in row] for row in U]


def reduce_binary(H):
    K = identity(2)
    while True:
        A = multiply(multiply(transpose(K), H), K)
        a, b, c = A[0][0], A[0][1], A[1][1]
        q = (2*b + a) // (2*a)
        if q:
            K = multiply(K, [[1, -q], [0, 1]])
            A = multiply(multiply(transpose(K), H), K)
            a, b, c = A[0][0], A[0][1], A[1][1]
        if c < a:
            K = multiply(K, [[0, 1], [1, 0]])
            continue
        if b < 0:
            K = multiply(K, [[1, 0], [0, -1]])
        return K


def split_basis(G, v):
    gv, d = matvec(G, v), dot(v, matvec(G, v))
    require(d > 0 and all(x % d == 0 for x in gv), "v does not split integrally")
    lam = [x//d for x in gv]
    U = complete_primitive(v)
    columns = transpose(U)
    for i in (1, 2):
        k = dot(lam, columns[i])
        columns[i] = [x-k*y for x, y in zip(columns[i], v)]
    B = transpose(columns)
    gram = multiply(multiply(transpose(B), G), B)
    H2 = [row[1:] for row in gram[1:]]
    K2 = reduce_binary(H2)
    K = [[1, 0, 0], [0, K2[0][0], K2[0][1]], [0, K2[1][0], K2[1][1]]]
    B = multiply(B, K)
    gram = multiply(multiply(transpose(B), G), B)
    return {"v": v, "norm": d, "projection_row": lam, "basis": B, "block_gram": gram,
            "type": "three_lines" if gram[1][2] == 0 else "line_and_indecomposable_binary"}


def decomposition(G, max_nodes=1_000_000):
    admit(G)
    bound = max(G[i][i] for i in range(3))
    vectors = sphere(G, bound, max_nodes)
    for v in sorted(vectors):
        if gcd(gcd(v[0], v[1]), v[2]) != 1:
            continue
        d = dot(v, matvec(G, v))
        if d and all(x % d == 0 for x in matvec(G, v)):
            return {"verdict": "decomposable", "search_bound": bound, "sphere_complete": True,
                    "split": split_basis(G, list(v))}
    return {"verdict": "indecomposable", "search_bound": bound, "sphere_complete": True}


def obtuse_superbase(G, max_steps=1_000_000):
    admit(G)
    integer(max_steps, "max_steps")
    require(max_steps >= 0, "negative Selling step limit")
    S = [[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    steps = []
    while True:
        pair = next(((i, j) for i, j in combinations(range(4), 2)
                     if dot(S[i], matvec(G, S[j])) > 0), None)
        if pair is None:
            break
        if len(steps) >= max_steps:
            raise EnumerationLimit("Selling reduction exceeded step limit")
        i, j = pair
        old = S[i][:]
        S[i] = [-x for x in old]
        for k in range(4):
            if k not in (i, j):
                S[k] = [x+y for x, y in zip(S[k], old)]
        steps.append([i, j])
    conorms = [[0 if i == j else -dot(S[i], matvec(G, S[j])) for j in range(4)] for i in range(4)]
    return {"vectors": S, "steps": steps, "conorms": conorms}


def connected(W, excluded=()):
    vertices = set(range(4)) - set(excluded)
    if not vertices:
        return True
    seen, pending = {min(vertices)}, [min(vertices)]
    while pending:
        i = pending.pop()
        for j in vertices - seen:
            if W[i][j] > 0:
                seen.add(j)
                pending.append(j)
    return seen == vertices


def graph_type(W):
    require(connected(W), "disconnected conorm graph")
    if any(not connected(W, (i,)) for i in range(4)):
        return "articulation"
    edges = sum(W[i][j] > 0 for i, j in combinations(range(4), 2))
    require(edges in (4, 5, 6), "unexpected graph without articulation")
    return {4: "cycle", 5: "diamond", 6: "complete"}[edges]


def shell_claim(G, a, b, superbase):
    """Symbolic complete shell from the four graph cases in RANK3-STRUCTURE.md."""
    require(dot(a, b) % 2 == 0, "even pair required")
    S, W = superbase["vectors"], superbase["conorms"]
    require(graph_type(W) != "articulation", "positive theorem needs no articulation")
    U = transpose(S[1:])
    ap = [int(x) % 2 for x in matvec(inverse(U), a)]
    A = [0] + ap
    odd = [i for i in range(4) if A[i]]
    if len(odd) > 2:
        A = [1-x for x in A]
        odd = [i for i in range(4) if A[i]]
    even = [i for i in range(4) if not A[i]]
    m = sum(W[i][j] for i in odd for j in even)
    case, xvectors = "constant", [[0]*4]
    if len(odd) == 1:
        case, xvectors = "singleton", []
        for sign in (-1, 1):
            x = [0]*4
            x[odd[0]] = sign
            xvectors.append(x)
    elif len(odd) == 2:
        wi, we = W[odd[0]][odd[1]], W[even[0]][even[1]]
        if wi and we:
            case, xvectors = "two_connected", []
            for sign in (-1, 1):
                x = [0]*4
                for i in odd:
                    x[i] = sign
                xvectors.append(x)
        else:
            if wi:
                odd, even = even, odd
                wi, we = we, wi
            case, xvectors = "diamond_minimum" if we else "cycle", []
            for signs in product((-1, 1), repeat=2):
                x = [0]*4
                for i, sign in zip(odd, signs):
                    x[i] = sign
                xvectors.append(x)
            if not we:
                for sign in (-1, 1):
                    x = [0]*4
                    x[even[1]] = 2*sign
                    for i in odd:
                        x[i] = sign
                    xvectors.append(x)
            elif all(dot(S[i], b) % 2 for i in odd):
                case, m, xvectors = "diamond_later", m+4*we, []
                for sign in (-1, 1):
                    x = [0]*4
                    x[even[1]] = 2*sign
                    for i in odd:
                        x[i] = sign
                    xvectors.append(x)
    vectors = [[sum(xi*vi[j] for xi, vi in zip(x, S)) for j in range(3)] for x in xvectors]
    require(all(all((yi-ai) % 2 == 0 for yi, ai in zip(y, a)) for y in vectors), "symbolic parity mismatch")
    require(all(dot(y, matvec(G, y)) == m for y in vectors), "symbolic norm mismatch")
    coeff = sum((-1)**(dot([(yi-ai)//2 for yi, ai in zip(y, a)], b) % 2) for y in vectors)
    require(coeff != 0, "symbolic coefficient unexpectedly zero")
    return {"norm4": m, "coefficient": coeff, "case": case,
            "selected_vectors": vectors, "complete_by": "rank3-obtuse-graph-shell-v1"}
