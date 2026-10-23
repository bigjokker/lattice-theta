"""Check saved cone maps against a fresh transcription of arXiv:2407.03513 p.5."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'bot_work/bot/coverage'))
import replay as R

OUT = Path(__file__).resolve().parent
vectors = [[int(i == j) for i in range(4)] for j in range(4)]
vectors += [[int(k == i)-int(k == j) for k in range(4)]
            for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]
rays = [[[a*b for b in v] for a in v] for v in vectors]
rays += [[[2,1,-1,-1],[1,2,-1,-1],[-1,-1,2,0],[-1,-1,0,2]]]
v = [1,1,-1,-1]
rays += [[[a*b for b in v] for a in v]]
assert R.is_pd(R.toQ(R.tovec(rays[10])))
pub = {'Delta':list(range(10)), 'Delta_prime':[0,1,2,3,5,6,7,8,9,10],
       'Delta_double_prime':[0,1,2,3,5,6,7,8,10,11]}
base = ROOT/'results/rank4-secondary-fan'
match = json.loads((base/'match-published.json').read_text())
rows = []
for name in ('fan.json','fan-exact.json'):
    fan = json.loads((base/name).read_text())
    assert len(fan['cones']) == len(match[name]) == 3
    used = []
    for cone, record in zip(fan['cones'],match[name]):
        U, label = record['U'], record['matched']
        assert all(type(x) is int for row in U for x in row) and abs(R.det_int(U)) == 1
        assert len(cone['rays']) == 10 and len(R.rref(cone['rays'])[1]) == 10
        actual = {tuple(R.tovec(R.congr(U,R.toQ(ray)))) for ray in cone['rays']}
        expected = {tuple(R.tovec(rays[i])) for i in pub[label]}
        assert actual == expected
        used.append(label)
    assert sorted(used) == sorted(pub)
    rows.append({'file':name,'bijective_exact_ray_match':True})
summary = json.loads((ROOT/'results/rank4-six-cells/cells-v3.json').read_text())['summary']
for name, h in summary['source_sha256'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == h
assert hashlib.sha256((ROOT/summary['fan']).read_bytes()).hexdigest() == summary['fan_sha256']
result = {'published_cones':rows,'source_hashes_match_handoff':True,
          'paper':'https://arxiv.org/pdf/2407.03513','page':5,
          'ray_transcription':'R1..R10 outer products; R11 primitive half of printed matrix; R12 outer product (1,1,-1,-1)'}
(OUT/'published-cones-and-hashes.json').write_text(json.dumps(result,indent=1)+'\n')
print(json.dumps(result))
