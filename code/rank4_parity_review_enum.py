"""Reviewer's independent re-enumeration of the class lists of
RANK4-PARITY-CONJECTURE.md (Tasks B, C, E) with a different reduction box
(review record: RANK4-PARITY-REVIEW.md).  Uses only rank4_feasibility and
exact_theta; nothing from rank4_parity_search.py or verify_rank4_parity.py.

Box.  Every positive definite rank-four lattice has a Minkowski-reduced basis:
G_kk = min Q(x) over integral x with gcd(x_k, ..., x_4) = 1, k = 1..4.  Hence
  (i)   G11 <= G22 <= G33 <= G44;
  (ii)  Q(x) >= G_kk for every x in {0, +-1}^4 with x_k = 1 and x_l = 0 for
        l > k (such x has gcd(x_k, ..., x_4) = 1); this contains |2 G_ij| <= G_ii;
  (iii) G11 G22 G33 G44 <= lambda_4 det G with lambda_4 = 4: in rank <= 4 a
        Minkowski-reduced basis attains the successive minima, and Minkowski's
        second theorem gives their product <= gamma_4^4 det = 4 det.  The
        constant is the --lam option so that a looser value serves as a control.
Signs: replacing e_j by -e_j (j = 2, 3, 4) flips G_1j, so G12, G13, G14 >= 0
may be assumed.  Candidates are bases; exact isometry
(rank4_feasibility.isometries, complete basis-image search) deduplicates,
after a cheap invariant (determinant and vector counts of all norms <= 40).
Decisions: a zero is an epsilon-one stabilizer element (exact_theta.
symmetry_data), a nonzero is a nonzero signed coefficient of a complete shell
(exact_theta.complete_shells) on the ladder norm4 100, 400, 1600; anything
else is reported as undecided.

Examples (from the workspace root; the output path must be fresh):
    python code/rank4_parity_review_enum.py --universe B --dmax 201 --output reports/rank4-parity/review-independent-enum-B201.txt
    python code/rank4_parity_review_enum.py --universe C --dmax 64  --output ...
    python code/rank4_parity_review_enum.py --universe E --dmin 25 --dmax 168 --output ...
"""
import argparse, sys, os, time, json
from itertools import product
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank4_feasibility import isometries, determinant
from exact_theta import complete_shells, symmetry_data, multiply, transpose, matvec, dot, inverse, identity, EnumerationLimit

BITS = [list(v) for v in product((0, 1), repeat=4)]
EVEN = [(a, b) for a in BITS if any(a) for b in BITS if any(b) and dot(a, b) % 2 == 0]
I4 = identity(4)
MI4 = [[-x for x in r] for r in I4]
INV_NORM = 40


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def det4(M):
    total = 0
    for j in range(4):
        minor = [[M[r][c] for c in range(4) if c != j] for r in range(1, 4)]
        total += (-1) ** j * M[0][j] * det3(minor)
    return total


def is_even(G):
    return all(G[i][i] % 2 == 0 for i in range(4))


def q(G, x):
    return dot(x, matvec(G, x))


K4 = [(x1, x2, x3, 1) for x1 in (-1, 0, 1) for x2 in (-1, 0, 1) for x3 in (-1, 0, 1)]
K3 = [(x1, x2, 1) for x1 in (-1, 0, 1) for x2 in (-1, 0, 1)]


def minkowski_candidates(dmin, dmax, even_only, odd_only, det_ok, lam):
    found = []
    bound = lam * dmax
    step = 2 if even_only else 1
    g11 = step
    while g11 ** 4 <= bound:
        g22 = g11
        while g11 * g22 ** 3 <= bound:
            g33 = g22
            while g11 * g22 * g33 ** 2 <= bound:
                g44 = g33
                while g11 * g22 * g33 * g44 <= bound:
                    for g12 in range(0, g11 // 2 + 1):
                        if g11 * g22 - g12 * g12 <= 0:
                            continue
                        for g13 in range(0, g11 // 2 + 1):
                            for g23 in range(-(g22 // 2), g22 // 2 + 1):
                                M3 = [[g11, g12, g13], [g12, g22, g23], [g13, g23, g33]]
                                if det3(M3) <= 0:
                                    continue
                                if any(g11 * x1 * x1 + g22 * x2 * x2 + g33 + 2 * (g12 * x1 * x2 + g13 * x1 + g23 * x2) < g33
                                       for x1, x2, _ in K3):
                                    continue
                                for g14 in range(0, g11 // 2 + 1):
                                    for g24 in range(-(g22 // 2), g22 // 2 + 1):
                                        for g34 in range(-(g33 // 2), g33 // 2 + 1):
                                            G = [[g11, g12, g13, g14], [g12, g22, g23, g24], [g13, g23, g33, g34], [g14, g24, g34, g44]]
                                            d = det4(G)
                                            if d < dmin or d > dmax or not det_ok(d):
                                                continue
                                            if odd_only and is_even(G):
                                                continue
                                            if any(q(G, list(x)) < g44 for x in K4):
                                                continue
                                            found.append(G)
                    g44 += step
                g33 += step
            g22 += step
        g11 += step
    return sorted(found)


def invariant(G, max_nodes=5_000_000):
    sh = complete_shells(G, [0] * 4, [0] * 4, 4 * INV_NORM, max_nodes)["shells"]
    return tuple((s["norm4"] // 4, s["vector_count"]) for s in sh if s["norm4"])


def order(T):
    P, k = T, 1
    while P != I4:
        P = multiply(P, T)
        k += 1
    return k


def splitting_involutions(group):
    return [T for T in group if T not in (I4, MI4) and multiply(T, T) == I4
            and all((x + y) % 2 == 0 for rT, rI in zip(T, I4) for x, y in zip(rT, rI))]


def eichler_components(G, max_nodes=20_000_000):
    R = max(G[i][i] for i in range(4))
    vecs = []
    complete_shells(G, [0] * 4, [0] * 4, 4 * R, max_nodes, lambda y, qq: vecs.append(tuple(v // 2 for v in y)))
    vecs = [v for v in vecs if any(v)]
    norm = {v: dot(v, matvec(G, v)) for v in vecs}
    by_norm = defaultdict(list)
    for v in vecs:
        by_norm[norm[v]].append(v)
    indec = []
    for v in vecs:
        if v < tuple(-x for x in v):
            continue
        Gv = matvec(G, v)
        if not any(dot(u, Gv) == n for n in by_norm if n < norm[v] for u in by_norm[n]):
            indec.append(v)
    parent = list(range(len(indec)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(len(indec)):
        Gi = matvec(G, indec[i])
        for j in range(i + 1, len(indec)):
            if dot(indec[j], Gi) != 0:
                parent[find(i)] = find(j)
    return len({find(i) for i in range(len(indec))})


def stabilizer_eps(G, a, b, group, duals):
    out = []
    for T, dual in zip(group, duals):
        if all((x - y) % 2 == 0 for x, y in zip(matvec(T, a), a)) and all((x - y) % 2 == 0 for x, y in zip(matvec(dual, b), b)):
            out.append((T, symmetry_data(G, a, b, T)["epsilon"]))
    return out


def decide(G, a, b, group, duals, ladder=(100, 400, 1600), max_nodes=20_000_000):
    st = stabilizer_eps(G, a, b, group, duals)
    wit = [T for T, e in st if e == 1]
    if wit:
        return "zero", sorted({order(T) for T in wit}), None
    for bound in ladder:
        try:
            sh = complete_shells(G, a, b, bound, max_nodes)["shells"]
        except EnumerationLimit:
            return "undecided", None, None
        nz = next((s for s in sh if s["signed_coefficient"]), None)
        if nz:
            return "nonzero", None, nz["norm4"]
    return "undecided", None, None


def singular_count(G):
    return sum(1 for x in BITS if any(x) and (q(G, x) // 2) % 2 == 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--universe", choices=("B", "C", "E"), required=True,
                    help="B: even lattices, odd det; C: even lattices, 4 | det; E: odd lattices, 4 does not divide det")
    ap.add_argument("--dmin", type=int, default=1)
    ap.add_argument("--dmax", type=int, required=True)
    ap.add_argument("--lam", type=int, default=4)
    ap.add_argument("--no-decide", action="store_true", help="class list only")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.exists(args.output):
        sys.exit("choose a fresh output path")
    out = open(args.output, "w", encoding="utf-8")

    def say(*s):
        print(*s, file=out, flush=True)

    def progress(msg):
        print(msg, file=sys.stderr, flush=True)
    even_only = args.universe in ("B", "C")
    odd_only = args.universe == "E"
    det_ok = {"B": lambda d: d % 2 == 1, "C": lambda d: d % 4 == 0, "E": lambda d: d % 4 != 0}[args.universe]
    label = {"B": "even lattices with odd determinant", "C": "even lattices with determinant divisible by 4",
             "E": "odd lattices with determinant not divisible by 4"}[args.universe]
    t0 = time.time()
    say(f"# Reviewer's independent re-enumeration: universe {args.universe} = isometry classes of positive definite integral "
        f"rank-four {label}, {args.dmin} <= det <= {args.dmax}; Minkowski box with lambda_4 = {args.lam}; exact-isometry dedup")
    cands = minkowski_candidates(args.dmin, args.dmax, even_only, odd_only, det_ok, args.lam)
    t_gen = time.time() - t0
    say(f"candidates (Minkowski-reduced Gram matrices in the box): {len(cands)}  generation {t_gen:.1f} s")
    progress(f"{len(cands)} candidates in {t_gen:.1f} s")
    reps, by_key, searches = [], defaultdict(list), 0
    for n, G in enumerate(cands):
        d = det4(G)
        key = (d, invariant(G))
        hit = None
        for i in by_key[key]:
            searches += 1
            if isometries(G, reps[i], 5_000_000):
                hit = i
                break
        if hit is None:
            by_key[key].append(len(reps))
            reps.append(G)
        if (n + 1) % 2000 == 0:
            progress(f"  candidates {n + 1}/{len(cands)} classes {len(reps)} searches {searches} {time.time() - t0:.0f}s")
    say(f"classes: {len(reps)}  (isometry searches {searches}; representatives of equal determinant and equal vector counts "
        f"through norm {INV_NORM} were compared by complete basis-image search)  {time.time() - t0:.1f} s")
    per_det = Counter(det4(G) for G in reps)
    say("classes per determinant: " + json.dumps(dict(sorted(per_det.items()))))
    if args.no_decide:
        out.close()
        return
    rows = []
    say("")
    say("| # | det | parity | |Aut| | sqrt(-I) | type (singular classes) | decomposable | P | zeros | witness orders | max first-nonzero norm4 | G |")
    say("|---:|---:|---|---:|---|---|---|---|---:|---|---:|---|")
    for i, G in enumerate(reps):
        d = det4(G)
        group = isometries(G, G, 5_000_000)
        duals = [transpose([[int(x) for x in r] for r in inverse(T)]) for T in group]
        inv = splitting_involutions(group)
        comps = eichler_components(G)
        if (len(inv) > 0) != (comps >= 2):
            say(f"DECOMPOSABILITY CRITERIA DISAGREE at class {i} det {d}: involutions {len(inv)} components {comps}")
        dec = len(inv) > 0
        sq = any(multiply(T, T) == MI4 for T in group)
        typ = "-"
        if is_even(G):
            sc = singular_count(G)
            typ = f"{'plus' if sc == 9 else 'minus' if sc == 5 else 'degenerate'} ({sc})"
        row = {"det": d, "even": is_even(G), "aut": len(group), "sqrt": sq, "decomposable": dec, "type": typ, "G": G}
        if not dec:
            verdicts, worders, firsts = Counter(), Counter(), []
            for a, b in EVEN:
                v, wo, first = decide(G, a, b, group, duals)
                verdicts[v] += 1
                if wo:
                    for o in wo:
                        worders[o] += 1
                if first is not None:
                    firsts.append(first)
            status = "holds" if verdicts["nonzero"] == 105 else "fails" if verdicts["zero"] else "undecided"
            row.update(P=status, zeros=verdicts["zero"], undecided=verdicts["undecided"], worders=dict(worders),
                       maxfirst=max(firsts) if firsts else None, firsts=firsts)
        else:
            row.update(P="fails (decomposable)", zeros=None, undecided=0, worders={}, maxfirst=None, firsts=[])
        rows.append(row)
        say(f"| {i} | {d} | {'even' if row['even'] else 'odd'} | {row['aut']} | {sq} | {typ} | {dec} | {row['P']} | "
            f"{row['zeros'] if row['zeros'] is not None else '-'} | {row['worders']} | {row['maxfirst']} | {G} |")
        if (i + 1) % 50 == 0:
            progress(f"  decided {i + 1}/{len(reps)} {time.time() - t0:.0f}s")
    indec = [r for r in rows if not r["decomposable"]]
    say("")
    say("## Summary")
    say(f"classes {len(rows)}; decomposable {len(rows) - len(indec)}; indecomposable {len(indec)}; "
        f"P holds {sum(1 for r in indec if r['P'] == 'holds')}; P fails {sum(1 for r in indec if r['P'] == 'fails')}; "
        f"undecided {sum(1 for r in indec if r['P'] == 'undecided')}")
    say("failing determinants (indecomposable): " + str(sorted(r["det"] for r in indec if r["P"] == "fails")))
    say("failing classes: " + json.dumps([{"det": r["det"], "aut": r["aut"], "zeros": r["zeros"], "witness_orders": r["worders"], "G": r["G"]}
                                           for r in indec if r["P"] == "fails"]))
    say("undecided characteristics: " + str(sum(r["undecided"] for r in indec)))
    say("indecomposable classes per determinant: " + json.dumps(dict(sorted(Counter(r["det"] for r in indec).items()))))
    say("decomposable classes per determinant: " + json.dumps(dict(sorted(Counter(r["det"] for r in rows if r["decomposable"]).items()))))
    say("|Aut| distribution of indecomposable classes: " + json.dumps(dict(sorted(Counter(r["aut"] for r in indec).items()))))
    say("type distribution of indecomposable classes: " + json.dumps(dict(Counter(r["type"].split(" ")[0] for r in indec))))
    say("type agrees with det mod 8 (plus iff 1, minus iff 5) on every even indecomposable class: "
        + str(all((r["type"].startswith("plus")) == (r["det"] % 8 == 1) and (r["type"].startswith("minus")) == (r["det"] % 8 == 5)
                  for r in indec if r["even"] and r["det"] % 2 == 1)))
    say("classes with an element squaring to -I: " + json.dumps([{"det": r["det"], "decomposable": r["decomposable"], "P": r["P"]} for r in rows if r["sqrt"]]))
    allfirst = Counter(f for r in indec for f in r["firsts"])
    say("first-nonzero norm4 distribution over all nonzero decisions of indecomposable classes: " + json.dumps(dict(sorted(allfirst.items()))))
    say(f"total nonzero decisions {sum(allfirst.values())}; zeros {sum(r['zeros'] or 0 for r in indec)}; "
        f"largest first-nonzero norm4 {max(allfirst) if allfirst else None}")
    say("zero witness orders over all zeros of indecomposable classes: " + json.dumps(dict(sorted(Counter(o for r in indec for o, n in r["worders"].items() for _ in range(n)).items()))))
    say(f"elapsed {time.time() - t0:.1f} s")
    out.close()


if __name__ == "__main__":
    main()
