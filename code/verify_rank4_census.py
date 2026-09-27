"""Independent boxes replay rank-four coverage, uniqueness, groups and coefficients."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells, dot,
                         integer, inverse, ldl, matrix, matvec, multiply, require,
                         symmetry_data, transpose, vector)
from verify_lattice_census import box_shells
from verify_rank4_feasibility import (BITS, box_isometries, cutoff_replay, det4,
                                      independent_candidates, modular_check)


def _replay(doc, max_nodes, backend):
    require(doc["schema"] == "work7-rank4-census-v1", "wrong census schema")
    require(backend in ("box", "ldl"), "wrong coefficient backend")
    H = integer(doc["domain"]["determinant_bound"], "determinant bound")
    require(integer(doc["domain"]["rank"], "rank") == 4 and 1 <= H <= 24, "wrong domain")
    bound = integer(doc["run"]["norm4_bound"], "norm4 bound")
    builder_nodes = integer(doc["run"]["max_nodes"], "builder nodes")
    integer(max_nodes, "replay nodes")
    require(bound >= 0 and builder_nodes > 0 and max_nodes > 0, "invalid limits")
    rows, generation = doc["lattices"], doc["generation"]
    for item in generation:
        matrix(item["G"], "candidate Gram", 4)
    require([g["G"] for g in generation] == independent_candidates(H), "candidate domain mismatch")
    for row in rows:
        G = matrix(row["G"], "representative Gram", 4)
        ldl(G)
        require(integer(row["determinant"], "determinant") == det4(G)
                and 1 <= det4(G) <= H, "wrong representative determinant")
    first_occurrences = []
    for item in generation:
        index = integer(item["representative"], "representative index")
        require(0 <= index < len(rows), "bad representative index")
        G, T = item["G"], matrix(item["T"], "generation transport", 4)
        require(abs(det4(T)) == 1 and multiply(multiply(transpose(T), rows[index]["G"]), T) == G,
                "invalid generation transport")
        if index not in first_occurrences:
            require(index == len(first_occurrences) and rows[index]["G"] == G,
                    "wrong representative retention order")
            first_occurrences.append(index)
    require(first_occurrences == list(range(len(rows))), "unused representative")
    for i, row in enumerate(rows):
        for previous in rows[:i]:
            if row["determinant"] == previous["determinant"]:
                require(not box_isometries(row["G"], previous["G"], max_nodes),
                        "duplicate isometry class")
    even, reports = Counter(), []
    for index, row in enumerate(rows):
        G = row["G"]
        group = box_isometries(G, G, max_nodes)
        for T in row["automorphisms"]:
            matrix(T, "stored automorphism", 4)
        require(row["automorphisms"] == group, "incomplete or incorrect full group")
        duals = [transpose(inverse(T)) for T in group]
        # Independent direct integral actions; cache only matrix-vector products.
        primal = [[matvec(T, a) for a in BITS] for T in group]
        dual = [[matvec(T, b) for b in BITS] for T in duals]
        require([m["b"] for m in row["modular_kernels"]] == [b for b in BITS if any(b)],
                "incomplete modular kernel coverage")
        for item in row["modular_kernels"]:
            modular_check(G, [0]*4, item["b"], item)
        require([(p["a"], p["b"]) for p in row["pairs"]] == [(a,b) for a in BITS for b in BITS],
                "characteristic domain mismatch")
        class_counts = Counter()
        for p in row["pairs"]:
            a, b = vector(p["a"], "a", 4), vector(p["b"], "b", 4)
            ai, bi = BITS.index(a), BITS.index(b)
            indices, signs = [], []
            for i in range(len(group)):
                if any((v-w) % 2 for v,w in zip(primal[i][ai], a)):
                    continue
                transported = dual[i][bi]
                if any((v-w) % 2 for v,w in zip(transported, b)):
                    continue
                shift = [(v-w)/2 for v,w in zip(transported,b)]
                require(all(v.denominator == 1 for v in shift), "nonintegral shift")
                indices.append(i)
                signs.append(int(dot(a,shift)) % 2)
            for value in p["stabilizer"]["indices"] + p["stabilizer"]["epsilon"]:
                integer(value, "stabilizer entry")
            require(p["stabilizer"] == {"indices": indices, "epsilon": signs},
                    "incomplete stabilizer or incorrect epsilon")
            verdict = p["verdict"]
            if verdict == "proved_zero":
                require(p["evidence"]["kind"] == "symmetry" and 1 in signs
                        and symmetry_data(G,a,b,p["evidence"]["T"])["epsilon"] == 1,
                        "zero lacks symmetry witness")
            else:
                require(1 not in signs, "unexplained verdict despite zero symmetry")
                if "enumeration" in p:
                    stored = p["enumeration"]
                    require(stored["complete"] is True
                            and integer(stored["norm4_bound"], "enumeration bound") == bound,
                            "wrong coefficient completeness")
                    shells = (box_shells(G,a,b,bound,max_nodes) if backend == "box"
                              else complete_shells(G,a,b,bound,max_nodes)["shells"])
                    require(shells == stored["shells"], "coefficient histogram mismatch")
                else:
                    require(p.get("enumeration_complete") is False and isinstance(p.get("reason"),str),
                            "missing coefficient evidence")
                    shells = None
                first = next((s for s in shells or [] if s["signed_coefficient"]), None)
                if "cutoff_certificate" in p:
                    cert = p["cutoff_certificate"]
                    require(cert["G"] == G and cert["a"] == a and cert["b"] == b,
                            "cutoff attached to wrong characteristic")
                    checked = cutoff_replay(cert,max_nodes)
                    require(checked["certificate_replay_complete"], "cutoff replay interrupted")
                    outcome = checked["verdict"]
                else:
                    outcome = "unresolved"
                if verdict == "proved_nonzero":
                    claim = p["evidence"]
                    require(first is not None and claim == {"kind":"shell", "norm4":first["norm4"],
                            "coefficient":first["signed_coefficient"]}, "nonzero lacks complete shell")
                elif verdict == "proved_zero_without_symmetry":
                    require(first is None and outcome == "proved_zero_by_modular_cutoff",
                            "false converse counterexample")
                elif verdict == "proved_nonzero_by_modular_replay":
                    require(first is None and outcome == "proved_nonzero", "invalid cutoff nonzero")
                else:
                    require(verdict == "unresolved" and first is None and outcome == "unresolved",
                            "wrong unresolved verdict")
            if dot(a,b) % 2 == 0:
                even[verdict] += 1
                class_counts[verdict] += 1
            else:
                require(verdict == "proved_zero", "odd pair must have symmetry zero")
        reports.append({"representative":index,"determinant":row["determinant"],
                        "group_order":len(group),"even_verdicts":dict(class_counts)})
    expected = {"candidate_bases":len(generation), "lattice_classes":len(rows),
                "even_pairs":136*len(rows), "odd_pairs":120*len(rows), "even_verdicts":dict(even)}
    for key in ("candidate_bases", "lattice_classes", "even_pairs", "odd_pairs"):
        integer(doc["summary"][key], key)
    for value in doc["summary"]["even_verdicts"].values():
        integer(value, "verdict count")
    require(doc["summary"] == expected, "wrong census summary")
    return {"verdict":"verified", "summary":expected, "lattices":reports,
            "domain_coverage_verified":True, "representative_uniqueness_verified":True,
            "full_groups_and_stabilizers_verified":True,
            "coefficient_backend":backend, "unresolved_pairs":even["unresolved"],
            "symmetry_converse_in_domain":False if even["proved_zero_without_symmetry"] else
                None if even["unresolved"] else True}


def replay(doc, max_nodes=1_000_000, backend="box"):
    try:
        return _replay(doc,max_nodes,backend)
    except EnumerationLimit as error:
        return {"verdict":"unresolved","reason":str(error),"symmetry_converse_in_domain":None}
    except (KeyError, TypeError, IndexError, ValueError, ZeroDivisionError) as error:
        raise InvalidCertificate("invalid rank-four census: "+str(error)) from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--backend",choices=("box","ldl"),default="box")
    parser.add_argument("--max-nodes",type=int,default=1_000_000)
    args = parser.parse_args()
    if args.output:
        require(not args.output.exists(), "choose a fresh replay path")
    raw = args.certificate.read_bytes()
    start = time.perf_counter()
    result = replay(json.loads(raw),args.max_nodes,args.backend)
    result.update(input_sha256=hashlib.sha256(raw).hexdigest(),seconds=time.perf_counter()-start)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({key:value for key,value in result.items() if key != "lattices"}))
    return 0 if result["verdict"] == "verified" else 2


if __name__ == "__main__":
    raise SystemExit(main())
