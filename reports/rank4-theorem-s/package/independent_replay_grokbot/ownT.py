"""Second pass: for every certified case, find T by Grok Bot's own complete isometry search (saved T NOT used), check D = Z E."""
import gzip, json, random, sys, time
from fractions import Fraction as F
import replay as r
fan = json.load(open(sys.argv[1])); cones = [c["rays"] for c in fan["cones"]]
rng = random.Random(7); t0 = time.time(); n_ok = 0; fail = []; tot = 0
for line in gzip.open(sys.argv[2], "rt"):
    rec = json.loads(line)
    if rec["d"] != "certified": continue
    tot += 1
    R = [cones[rec["cone"]][j] for j in rec["subset"]]; n = len(R); E = rec["E"]; tau = tuple(rec["tau"])
    mu = [F(x) for x in rec["mu"]]; Nb = r.nullspace(E, n)
    while True:
        g = [F(rng.randint(-1000, 1000), 1000) for _ in Nb]
        mg = [mu[j] + F(1, 4)*sum(g[i]*Nb[i][j] for i in range(len(Nb))) for j in range(n)]
        if all(x > 0 for x in mg): break
    Hg = r.toQ([sum(mg[j]*R[j][t] for j in range(n)) for t in range(10)])
    found = None
    for T in r.isometries(Hg, tau):
        if r.mod2_perm(T) != tau: continue
        D = [[F(0)]*n for _ in range(10)]
        for j in range(n):
            dv = r.tovec(r.congr(T, r.toQ(R[j])))
            for t in range(10): D[t][j] = dv[t] - F(R[j][t])
        if r.solve_left(D, E) is not None:
            found = T; break
    if found: n_ok += 1
    else: fail.append([rec["cone"], rec["subset"], rec["tau"]])
print(json.dumps({"certified_cases": tot, "own_T_with_DZE": n_ok, "failures": fail, "sec": round(time.time()-t0, 1)}))
