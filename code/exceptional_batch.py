"""Batch exact shell coefficients and compact symmetry certificates for E7/E8."""

from collections import Counter, defaultdict, deque
from fractions import Fraction
import hashlib
import json

from exact_theta import (complete_shells, identity, integer, inverse, ldl, matrix,
                         matvec, multiply, require, symmetry_data, transpose, verify)


def root_gram(rank):
    require(rank in (6, 7, 8), "exceptional rank must be 6, 7 or 8")
    G = [[2 * int(i == j) for j in range(rank)] for i in range(rank)]
    for i in range(rank - 2):
        G[i][i + 1] = G[i + 1][i] = -1
    G[2][-1] = G[-1][2] = -1
    return G


def bits(mask, rank):
    return [(mask >> i) & 1 for i in range(rank)]


def mask(vector):
    return sum((x % 2) << i for i, x in enumerate(vector))


def walsh(histogram, rank):
    """All sign sums in O(rank*2^rank) additions; no normalization."""
    values = [histogram.get(i, 0) for i in range(1 << rank)]
    step = 1
    while step < len(values):
        for base in range(0, len(values), 2 * step):
            for i in range(base, base + step):
                x, y = values[i], values[i + step]
                values[i], values[i + step] = x + y, x - y
        step *= 2
    return values


def ldl_histograms(G, bound, max_nodes=1_000_000):
    rank = len(G)
    histograms = defaultdict(Counter)
    node_count = vector_count = 0
    for am in range(1 << rank):
        a = bits(am, rank)
        private = defaultdict(Counter)

        def observe(y, N):
            xm = sum((((yi - ai) // 2) % 2) << i for i, (yi, ai) in enumerate(zip(y, a)))
            private[N][xm] += 1

        # A failure never publishes the partially accumulated private histogram.
        enumeration = complete_shells(G, a, [0] * rank, bound, max_nodes, observe)
        node_count += enumeration["visited_nodes"]
        for shell in enumeration["shells"]:
            N, count = shell["norm4"], shell["vector_count"]
            require(sum(private[N].values()) == count, "observer count mismatch")
            histograms[(am, N)] = private[N]
            vector_count += count
    return histograms, {"ldl_nodes": node_count, "lattice_vectors": vector_count,
                        "shifted_enumerations": 1 << rank, "norm4_bound": bound}


def histogram_hash(histograms):
    rows = [[am, N, sorted(histogram.items())] for (am, N), histogram in sorted(histograms.items())]
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode("ascii")).hexdigest()


def reflections(G):
    result = []
    for i in range(len(G)):
        S = identity(len(G))
        S[i] = [x - y for x, y in zip(S[i], G[i])]
        result.append(S)
    return result


def action_tables(G, generators):
    rank, size = len(G), 1 << len(G)
    tables = []
    for S in generators:
        matrix(S, "generator", rank)
        require(multiply(multiply(transpose(S), G), S) == G, "generator is not an isometry")
        # These proof words use self-inverse generators so path reversal is inversion.
        require(multiply(S, S) == identity(rank), "generator must be an involution")
        dual = transpose(S)
        primal = [mask(matvec(S, bits(am, rank))) for am in range(size)]
        dual_action, mu_masks = [], []
        for bm in range(size):
            raw = matvec(dual, bits(bm, rank))
            canonical = [x % 2 for x in raw]
            dual_action.append(mask(canonical))
            mu_masks.append(mask([(x - y) // 2 for x, y in zip(raw, canonical)]))
        tables.append((primal, dual_action, mu_masks))
    return tables


def step(tables, am, bm, gi):
    primal, dual, mu = tables[gi]
    aa, bb = primal[am], dual[bm]
    return aa, bb, (aa & mu[bm]).bit_count() % 2


def path_word(nodes, index):
    word = []
    while "parent" in nodes[index]:
        word.append(nodes[index]["generator"])
        index = nodes[index]["parent"]
    return list(reversed(word))


def word_matrix(generators, word, rank):
    T = identity(rank)
    for gi in word:
        integer(gi, "word generator")
        require(0 <= gi < len(generators), "word generator out of range")
        T = multiply(generators[gi], T)
    return T


def symmetry_forest(G, candidates):
    """Negative signed cycles produce root witnesses; tree paths transport zeros."""
    rank = len(G)
    generators = reflections(G)
    tables = action_tables(G, generators)
    candidates = set(candidates)
    nodes, seeds, indices, signs = [], [], {}, []
    for root in sorted(candidates):
        if root in indices:
            continue
        root_index = len(nodes)
        indices[root] = root_index
        nodes.append({"a_mask": root[0], "b_mask": root[1]})
        signs.append(0)
        queue = deque([root_index])
        negative_word = None
        while queue:
            source = queue.popleft()
            am, bm = nodes[source]["a_mask"], nodes[source]["b_mask"]
            for gi in range(rank):
                aa, bb, edge_sign = step(tables, am, bm, gi)
                require((aa, bb) in candidates, "zero-coefficient candidates are not invariant")
                target = indices.get((aa, bb))
                predicted = signs[source] ^ edge_sign
                if target is None:
                    target = len(nodes)
                    indices[(aa, bb)] = target
                    nodes.append({"a_mask": aa, "b_mask": bb, "parent": source, "generator": gi})
                    signs.append(predicted)
                    queue.append(target)
                elif negative_word is None and predicted != signs[target]:
                    negative_word = path_word(nodes, source) + [gi] + list(reversed(path_word(nodes, target)))
        require(negative_word is not None, "candidate component has no symmetry proof; unresolved")
        T = word_matrix(generators, negative_word, rank)
        require(symmetry_data(G, bits(root[0], rank), bits(root[1], rank), T)["epsilon"] == 1,
                "negative cycle did not produce epsilon=1")
        nodes[root_index]["seed"] = len(seeds)
        seeds.append({"a_mask": root[0], "b_mask": root[1], "T": T, "word": negative_word})
    return generators, seeds, nodes, indices


def validate_pack(pack):
    require(isinstance(pack, dict) and pack.get("schema") == "work7-theta-classification-v1", "unknown schema")
    G = matrix(pack.get("G"), "G")
    rank, size = len(G), 1 << len(G)
    require(rank in (7, 8) and G == root_gram(rank), "specified E7/E8 Gram required")
    _, pivots = ldl(G)
    determinant = Fraction(1)
    for pivot in pivots:
        determinant *= pivot
    require(pack.get("representatives") == "binary-primal-and-dual-masks", "unsupported representatives")
    bound = integer(pack.get("norm4_bound"), "norm4_bound")
    require(bound >= 0, "bound must be nonnegative")

    def key(obj):
        require(isinstance(obj, dict), "record must be an object")
        am, bm = integer(obj.get("a_mask"), "a_mask"), integer(obj.get("b_mask"), "b_mask")
        require(0 <= am < size and 0 <= bm < size, "characteristic mask out of range")
        require((am & bm).bit_count() % 2 == 0, "odd class in even classification")
        return am, bm

    generators = pack.get("generators")
    require(isinstance(generators, list) and generators, "generators required")
    tables = action_tables(G, generators)
    seeds = pack.get("seeds")
    require(isinstance(seeds, list), "seeds must be a list")
    for seed in seeds:
        am, bm = key(seed)
        require(isinstance(seed.get("word"), list), "seed word must be a list")
        require(word_matrix(generators, seed["word"], rank) == seed.get("T"), "seed word does not equal T")
        verify({"schema": "work7-theta-v1", "G": G, "a": bits(am, rank), "b": bits(bm, rank),
                "evidence": {"kind": "symmetry", "T": seed["T"]}})
    nodes = pack.get("zero_nodes")
    require(isinstance(nodes, list), "zero_nodes must be a list")
    node_keys, used_seeds = set(), set()
    for i, node in enumerate(nodes):
        am, bm = key(node)
        require((am, bm) not in node_keys, "duplicate symmetry node")
        node_keys.add((am, bm))
        if "parent" in node:
            require("seed" not in node, "node cannot have both parent and seed")
            parent, gi = integer(node["parent"], "parent"), integer(node.get("generator"), "generator")
            require(0 <= parent < i, "parent must refer to an earlier node")
            require(0 <= gi < len(generators), "generator index out of range")
            previous = nodes[parent]
            aa, bb, _ = step(tables, previous["a_mask"], previous["b_mask"], gi)
            require((aa, bb) == (am, bm), "symmetry transport reaches the wrong characteristic")
        else:
            si = integer(node.get("seed"), "seed")
            require(0 <= si < len(seeds), "seed index out of range")
            require(key(seeds[si]) == (am, bm), "root does not match its symmetry seed")
            used_seeds.add(si)
    require(used_seeds == set(range(len(seeds))), "unused seed")
    classes = pack.get("classes")
    require(isinstance(classes, list), "classes must be a list")
    expected = {(am, bm) for am in range(size) for bm in range(size) if (am & bm).bit_count() % 2 == 0}
    seen, referenced_nodes = set(), set()
    for entry in classes:
        am, bm = key(entry)
        require((am, bm) not in seen, "duplicate class")
        seen.add((am, bm))
        evidence = entry.get("evidence")
        require(isinstance(evidence, dict), "evidence required")
        if evidence.get("kind") == "symmetry":
            ni = integer(evidence.get("node"), "node")
            require(0 <= ni < len(nodes), "node index out of range")
            require(key(nodes[ni]) == (am, bm), "class references the wrong symmetry node")
            referenced_nodes.add(ni)
        elif evidence.get("kind") == "shell":
            N = integer(evidence.get("norm4"), "norm4")
            c = integer(evidence.get("coefficient"), "coefficient")
            require(0 <= N <= bound and c != 0, "nonzero shell within bound required")
        else:
            require(False, "unknown evidence kind")
    require(seen == expected, f"incomplete class coverage: expected {len(expected)}, found {len(seen)}")
    require(referenced_nodes == set(range(len(nodes))), "unused symmetry node")
    return G, rank, bound, int(determinant)


def structural_checks(pack):
    """Check the reported E7 glue and E8 quadratic-space characterizations."""
    G, rank = pack["G"], len(pack["G"])
    invG = inverse(G)
    if rank == 7:
        require(all((2 * x).denominator == 1 for row in invG for x in row), "2L* must lie in L")
        failures = 0
        for entry in pack["classes"]:
            if entry["evidence"]["kind"] == "symmetry":
                coords = matvec(invG, bits(entry["b_mask"], rank))
                require(any(x.denominator != 1 for x in coords), "E7 zero class has b in L")
                failures += 1
        return {"even_zeros_outside_primal_lattice": failures, "two_dual_lattice_is_in_primal": True}
    require(all(x.denominator == 1 for row in invG for x in row), "E8 must be unimodular")
    invG = [[int(x) for x in row] for row in invG]

    def q(am):
        v = bits(am, rank)
        norm = sum(x * y for x, y in zip(v, matvec(G, v)))
        require(norm % 2 == 0, "E8 must be even")
        return (norm // 2) % 2

    planes = set()
    for entry in pack["classes"]:
        am = entry["a_mask"]
        cm = mask(matvec(invG, bits(entry["b_mask"], rank)))
        singular_plane = am != 0 and cm != 0 and am != cm and q(am) == q(cm) == 0
        require(singular_plane == (entry["evidence"]["kind"] == "symmetry"),
                "E8 classification disagrees with the singular-plane criterion")
        if singular_plane:
            require(q(am ^ cm) == 0, "plane is not totally singular")
            planes.add(tuple(sorted((am, cm, am ^ cm))))
    # Exhibit a totally singular four-space to certify plus type, rather than
    # inferring its type solely from the singular-vector count.
    singular_basis, span = [], {0}
    for am in range(1, 1 << rank):
        if am in span or q(am):
            continue
        if any((am & mask(matvec(G, bits(v, rank)))).bit_count() % 2 for v in singular_basis):
            continue
        singular_basis.append(am)
        span |= {am ^ v for v in list(span)}
        if len(singular_basis) == 4:
            break
    require(len(span) == 16 and all(q(v) == 0 for v in span), "no totally singular four-space found")
    return {"totally_singular_planes": len(planes), "singular_vectors_including_zero":
            sum(q(am) == 0 for am in range(1 << rank)), "quadratic_criterion_verified": True,
            "totally_singular_four_space_basis_masks": singular_basis,
            "one_zero_orbit_under_supplied_generators": len(pack["seeds"]) == 1}


def export_certificate(pack, entry):
    """Expand a forest certificate to an ordinary G,a,b,T certificate."""
    rank = len(pack["G"])
    evidence = entry["evidence"]
    if evidence["kind"] == "symmetry":
        ni = evidence["node"]
        word = path_word(pack["zero_nodes"], ni)
        root = ni
        while "parent" in pack["zero_nodes"][root]:
            root = pack["zero_nodes"][root]["parent"]
        seed = pack["seeds"][pack["zero_nodes"][root]["seed"]]
        U = word_matrix(pack["generators"], word, rank)
        Ui = [[int(x) for x in row] for row in inverse(U)]
        evidence = {"kind": "symmetry", "T": multiply(multiply(U, seed["T"]), Ui)}
    return {"schema": "work7-theta-v1", "name": f"E{rank}: a={entry['a_mask']}, b={entry['b_mask']}",
            "G": pack["G"], "a": bits(entry["a_mask"], rank), "b": bits(entry["b_mask"], rank),
            "evidence": evidence}
