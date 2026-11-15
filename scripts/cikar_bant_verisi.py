# -*- coding: utf-8 -*-
"""
Extract per-registry band data from the registry-scan outputs
(gorev-a worktree 02_registry/runs_*/<mat>_<reg>_bnd.out) into
../data/registry_bantlar.json.

For each registry: eigenvalues at Gamma, K, K-0.02, K-0.04 (crystal, along
K->Gamma) for bands 1-2 (metal-s semicore, one per layer), VBM (36) and CBM
(37), plus the registry shift. This is the only step that touches files
outside preprint/; everything downstream reads data/ only.
"""
import json, os, re, sys

SRC = os.environ.get("REGISTRY_DIR", "02_registry")  # directory holding runs_<mat>/ and registry_tarama_<mat>.csv
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "registry_bantlar.json")
REGS = ["MMp", "d13_23", "d23_13"] + [f"ws{i}" for i in range(9)]
IV = 36   # 72 valence electrons (Se 3d semicore included) -> VBM = band 36


def blocks(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    if "JOB DONE" not in txt:
        raise RuntimeError(f"unfinished: {path}")
    parts = re.split(r"bands \(ev\):", txt)[1:]
    return [[float(x) for x in re.findall(r"-?\d+\.\d+", p.split("occupation")[0])]
            for p in parts]


def main():
    alat = {}
    data = {}
    for mat in ("WSSe", "MoSSe"):
        csv = os.path.join(SRC, f"registry_tarama_{mat}.csv")
        shifts = {}
        for line in open(csv).read().strip().splitlines()[1:]:
            c = line.split(",")
            shifts[c[0]] = (float(c[1]), float(c[2]))
        data[mat] = {}
        for reg in REGS:
            p = os.path.join(SRC, f"runs_{mat}", f"{mat}_{reg}_bnd.out")
            b = blocks(p)
            assert len(b) == 5, (p, len(b))
            G, K, M, Kd1, Kd2 = b
            if mat not in alat:
                m = re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)",
                              open(p, encoding="utf-8", errors="replace").read())
                alat[mat] = float(m.group(1)) * 0.529177210903  # bohr -> A
            pick = lambda e: dict(sc1=e[0], sc2=e[1], vb=e[IV - 1], cb=e[IV])
            data[mat][reg] = dict(d=shifts[reg], G=pick(G), K=pick(K),
                                  Kd1=pick(Kd1), Kd2=pick(Kd2))
    out = dict(aciklama=__doc__.strip().splitlines()[0],
               kaynak="gorev-a-janus-d798be/kol_A_janus/02_registry/runs_*",
               a0_ang=alat, dk_Kd1_2pi_over_a=0.04, dk_Kd2_2pi_over_a=0.08,
               bant_vbm=IV, veri=data)
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", os.path.normpath(OUT), {m: len(v) for m, v in data.items()})


if __name__ == "__main__":
    main()
