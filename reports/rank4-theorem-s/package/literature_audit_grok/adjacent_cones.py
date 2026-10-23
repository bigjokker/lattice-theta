"""Cross two facets of Delta = cone{R1..R10} and name the adjacent rays."""
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import rank4_secondary_fan as SF

PUB = {
    "R1": [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    "R2": [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    "R3": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
    "R4": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]],
    "R5": [[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    "R6": [[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]],
    "R7": [[1, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 1]],
    "R8": [[0, 0, 0, 0], [0, 1, -1, 0], [0, -1, 1, 0], [0, 0, 0, 0]],
    "R9": [[0, 0, 0, 0], [0, 1, 0, -1], [0, 0, 0, 0], [0, -1, 0, 1]],
    "R10": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, -1], [0, 0, -1, 1]],
    "R11": [[4, 2, -2, -2], [2, 4, -2, -2], [-2, -2, 4, 0], [-2, -2, 0, 4]],
    "R12": [[1, 1, -1, -1], [1, 1, -1, -1], [-1, -1, 1, 1], [-1, -1, 1, 1]],
}
QSET = {tuple(SF.to_q(Q)): name for name, Q in PUB.items()}


def name_ray(q):
    k = tuple(SF.primitive(q))
    return QSET.get(k, k)


def cross(drop):
    names = [f"R{i}" for i in range(1, 11) if f"R{i}" != drop]
    fpt = [0] * 10
    for n in names:
        for t, v in enumerate(SF.to_q(PUB[n])):
            fpt[t] += v
    # central of Delta is sum of R1..R10, already primitive
    cq = [0] * 10
    for i in range(1, 11):
        for t, v in enumerate(SF.to_q(PUB[f"R{i}"])):
            cq[t] += v
    t = F(1, 8)
    Qn = SF.to_Q(SF.primitive([x + t * (x - y) for x, y in zip(fpt, cq)]))
    print("drop", drop, "generic", Qn, "pd", SF.positive_definite(Qn), flush=True)
    cone = SF.secondary_cone(Qn)
    rays = [name_ray(r) for r in cone["rays"]]
    print("n", len(cone["rays"]), "central_det", int(SF.det(SF.central(cone))), flush=True)
    for r in rays:
        print(" ", r, flush=True)


if __name__ == "__main__":
    cross("R5")
    cross("R1")
