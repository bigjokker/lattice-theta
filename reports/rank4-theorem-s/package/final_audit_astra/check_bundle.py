"""Direct audit of every stored Grok Bot disposition, independent of its LP solver.

Uses Grok Bot's reviewed exact short-vector enumeration for fresh minima;
does not call Claude's classifier, cdd, simplex, or solve_left.
"""
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'bot_work/bot/coverage'))
import replay as R

OUT = Path(__file__).resolve().parent


def face_data(task):
    cid, mask, rays = task
    ids = tuple(j for j in range(10) if mask >> j & 1)
    rs = [rays[j] for j in ids]
    q = R.toQ([sum(r[t] for r in rs) for t in range(10)])
    if not R.is_pd(q):
        return (cid, ids), None
    sets, _ = R.minimal_sets(q)
    functions = {}
    for k, vectors in sets.items():
        values = {tuple(R.rayval(ray, v) for ray in rs) for v in vectors}
        assert len(values) == 1
        functions[k] = next(iter(values))
    return (cid, ids), (rs, sets, functions)


def main():
    start = time.perf_counter()
    fanpath = ROOT/'results/rank4-secondary-fan/fan.json'
    certpath = ROOT/'bot_work/bot/coverage/dispositions.jsonl.gz'
    cellpath = ROOT/'results/rank4-six-cells/cells.json'
    fan = json.loads(fanpath.read_text())
    assert len(fan['cones']) == 3 and all(len(c['rays']) == 10 for c in fan['cones'])
    tasks = [(cid, mask, cone['rays']) for cid, cone in enumerate(fan['cones']) for mask in range(1, 1024)]
    with get_context('spawn').Pool(16) as pool:
        faces = dict(pool.imap_unordered(face_data, tasks, chunksize=8))
    assert len(faces) == 3069
    faces = {k:v for k,v in faces.items() if v is not None}
    assert len(faces) == 2394
    print('Fresh exact minimal sets: 2394 positive definite faces', flush=True)
    taus = set(R.markings())
    assert len(taus) == 210
    seen = {k:set() for k in faces}
    counts = Counter()
    certified = set()
    with gzip.open(certpath, 'rt') as source:
        for line in source:
            rec = json.loads(line)
            key = (rec['cone'], tuple(rec['subset']))
            tau = tuple(rec['tau'])
            assert key in faces and tau in taus and tau not in seen[key]
            seen[key].add(tau)
            rays, mins, functions = faces[key]
            pairs = R.pairs_of(tau)
            assert len(pairs) == 6
            kind = rec['d']
            if kind == 'multiplicity_mismatch':
                a, b = rec['pair']
                assert tuple(sorted((a,b))) in pairs
                assert sorted(map(tuple,rec['Ma'])) == mins[a]
                assert sorted(map(tuple,rec['Mb'])) == mins[b]
                assert len(mins[a]) != len(mins[b])
            else:
                assert all(len(mins[a]) == len(mins[b]) for a,b in pairs)
                n = len(rays)
                E = [[functions[a][j]-functions[b][j] for j in range(n)] for a,b in pairs]
                mu = list(map(F, rec['mu']))
                assert len(mu) == n and all(v >= 0 for v in mu)
                assert all(sum(e[j]*mu[j] for j in range(n)) == 0 for e in E)
                if kind in ('lower_face','infeasible'):
                    y = list(map(F, rec['y']))
                    assert len(y) == len(E)
                    z = [sum(E[i][j]*y[i] for i in range(len(E))) for j in range(n)]
                    assert all(v >= 0 for v in z) and any(v > 0 for v in z)
                    assert z == list(map(F, rec['Ety']))
                    support = [rec['subset'][j] for j in range(n) if mu[j] > 0]
                    assert support == rec['support']
                    assert all(z[j] > 0 for j in range(n) if mu[j] == 0)
                    assert (kind == 'lower_face') == bool(support)
                else:
                    assert kind == 'certified' and all(v > 0 for v in mu)
                    assert E == rec['E']
                    T = rec['T']
                    assert len(T) == 4 and all(len(row) == 4 and all(type(v) is int for v in row) for row in T)
                    assert abs(R.det_int(T)) == 1 and R.mod2_perm(T) == tau
                    Z = [list(map(F,row)) for row in rec['Z']]
                    assert len(Z) == 10 and all(len(row) == len(E) for row in Z)
                    for j, ray in enumerate(rays):
                        moved = R.tovec(R.congr(T,R.toQ(ray)))
                        for t in range(10):
                            assert moved[t]-ray[t] == sum(Z[t][p]*E[p][j] for p in range(len(E)))
                    certified.add((key[0],key[1],tau))
            counts[kind] += 1
    assert all(value == taus for value in seen.values())
    saved = json.loads(cellpath.read_text())['cells']
    assert len(saved) == 5256
    assert certified == {(c['cone'],tuple(c['subset']),tuple(c['tau'])) for c in saved}
    assert counts == {'multiplicity_mismatch':477492,'lower_face':19992,'certified':5256}
    paths = [fanpath,certpath,cellpath,Path(R.__file__),Path(__file__)]
    result = {'faces':3069,'pd_faces':2394,'markings':210,'cases':sum(counts.values()),
              'counts':dict(counts),'complete_unique_coverage':True,'all_certificates_checked':True,
              'certified_set_matches':True,'wall_seconds':time.perf_counter()-start,
              'sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (OUT/'bundle-audit.json').write_text(json.dumps(result,indent=1)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'sha256'}),flush=True)


if __name__ == '__main__':
    main()
