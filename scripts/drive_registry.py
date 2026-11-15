"""Asama 2 — R-istif Janus cift katman registry taramasi (tMoSSe / tWSSe).

Her registry kaymasi delta icin:
  1. z-yalniz gevsetme (duzlem-ici kilitli, if_pos 0 0 1) — PBE+D3, dipfield
  2. ayni yogunlukla G+K+M+K-civari bant noktalari (SOC'siz) -> bant
     kenarlari, K etkin kutleleri, K'daki katmanlar-arasi yarilma
Cikti: registry_tarama_<mat>.csv + registry_ham_<mat>.json

Kullanim: python drive_registry.py [MoSSe|WSSe] [--probe]
  --probe: yalniz MM' noktasi (zamanlama testi)
"""
import os
import sys
import csv
import json
import re as _re
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from common import (hex_cell, bilayer_atoms, pw_input, run_pw, parse_pw,
                    parse_bands, RY_EV)

HERE = os.path.dirname(os.path.abspath(__file__))
MONO = os.path.join(HERE, "..", "01_monolayer")

PROD = json.load(open(os.path.join(MONO, "prod_params.json")))
ECUT = PROD["ecutwfc"]
# Cift katman gevsetme k=6: monolayer conv'da k6-k12 E farki 1 meV/hucre
# (registry FARKLARINDA buyuk olcude sistematik, iptal olur), dipol farki
# 0.0005 D. Kenar enerjileri bantlarda kesin k-noktalarinda hesaplanir.
KREL = PROD.get("kgrid_relax_bilayer", 6)
A0 = json.load(open(os.path.join(MONO, "a0_sonuc.json")))

C_BI = 30.0  # cift katman icin hucre yuksekligi (vakum > 20 A)
MATERIALS = {"MoSSe": ("Mo", "S", "Se"), "WSSe": ("W", "S", "Se")}

# Registry noktalari: yuksek simetri + 9 ic nokta (kristal koord)
HS = {"MMp": (0.0, 0.0), "d13_23": (1 / 3, 2 / 3), "d23_13": (2 / 3, 1 / 3)}
GRID9 = [(i / 3 + 1 / 6, j / 3 + 1 / 6) for i in range(3) for j in range(3)]
DELTAS = list(HS.items()) + [(f"ws{n}", d) for n, d in enumerate(GRID9)]


def reg_tag(mat, name):
    return f"{mat}_{name}"


def seed_zs(runs, mat, skip_name):
    """Bitmis bir registry kosusundan z-tohumlari (tur, z) listesi."""
    import glob as _g
    for fp in sorted(_g.glob(os.path.join(runs, f"{mat}_*_rlx.out")),
                     key=os.path.getmtime, reverse=True):
        if skip_name in os.path.basename(fp):
            continue
        r = parse_pw(fp)
        if r.get("ok") and r.get("final_coords_block"):
            zs = []
            for line in r["final_coords_block"].splitlines():
                p = line.split()
                if len(p) >= 4 and p[0] in ("Mo", "W", "S", "Se"):
                    zs.append(float(p[3]))
            if len(zs) == 6:
                return zs
    return None


def one_registry(mat, name, delta, runs):
    metal, xb, xt = MATERIALS[mat]
    a0 = A0[mat]["a0_ang"]
    tag = reg_tag(mat, name)
    res = {"registry": name, "delta": delta}

    # 1) z-gevsetme
    rtag = tag + "_rlx"
    outp = os.path.join(runs, rtag + ".out")
    if not (os.path.exists(outp) and parse_pw(outp).get("ok")):
        cell = hex_cell(a0, C_BI)
        atoms = bilayer_atoms(metal, xb, xt, C_BI, a0, delta_frac=delta)
        # onceki registry'nin z'leriyle tohumla (xy registry'ye ozgu kalir)
        zs = seed_zs(runs, mat, rtag)
        if zs:
            atoms = [(s, x, y, z) for (s, x, y, _), z in zip(atoms, zs)]
        txt = pw_input(f"reg_{tag}", cell, atoms, calculation="relax",
                       ecutwfc=ECUT, kgrid=(KREL, KREL, 1), soc=False,
                       dipfield=True, vdw="dft-d3", conv_thr=1e-8,
                       forc_conv_thr=1e-4, constrain_xy=True)
        txt = txt.replace("disk_io = 'low'", "disk_io = 'medium'")
        open(os.path.join(runs, rtag + ".in"), "w", newline="\n").write(txt)
        ok, _ = run_pw(runs, rtag + ".in", np=PROD.get("np", 8),
                       npool=PROD.get("npool", 4), timeout=4 * 3600,
                       keep_outdir=True)
        if not ok:
            print(f"[HATA] {rtag} gevsetme basarisiz", flush=True)
            return res
    r = parse_pw(outp)
    res.update(etot_ry=r.get("etot_ry"), dipole_debye=r.get("dipole_debye"),
               pot_amp_ry=r.get("pot_amp_ry"), wall_rlx=r.get("wall"))
    # gevsemis z'lerden katman araligi (arayuz kalkojen-kalkojen)
    if r.get("final_coords_block"):
        zs = []
        for line in r["final_coords_block"].splitlines():
            p = line.split()
            if len(p) >= 4 and p[0] in ("Mo", "W", "S", "Se"):
                zs.append((p[0], float(p[3])))
        if len(zs) == 6:
            res["z_crystal"] = zs
            res["d_interface_ang"] = (zs[4][1] - zs[2][1]) * C_BI

    # 2) bant noktalari: G, K, M (SOC'siz, yogunluk yeniden kullanilir)
    btag = tag + "_bnd"
    boutp = os.path.join(runs, btag + ".out")
    if not (os.path.exists(boutp) and "JOB DONE" in open(boutp, errors="ignore").read()):
        cell = hex_cell(a0, C_BI)
        # gevsemis koordinatlari dogrudan relax ciktisindan al
        pr = parse_pw(outp)
        atoms = []
        for line in pr["final_coords_block"].splitlines():
            p = line.split()
            if len(p) >= 4 and p[0] in ("Mo", "W", "S", "Se"):
                atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        # G, K, M + K civarinda 2 nokta (K->G yonu; parabolik kutle fiti)
        kexp = ("plain", [(0.0, 0.0, 0.0, 1.0),
                          (1 / 3, 1 / 3, 0.0, 1.0),
                          (0.5, 0.0, 0.0, 1.0),
                          (1 / 3 - 0.02, 1 / 3 - 0.02, 0.0, 1.0),
                          (1 / 3 - 0.04, 1 / 3 - 0.04, 0.0, 1.0)])
        txt = pw_input(f"reg_{tag}", cell, atoms, calculation="bands",
                       ecutwfc=ECUT, soc=False, dipfield=True, vdw="dft-d3",
                       nbnd=48, kpoints_explicit=kexp, conv_thr=1e-8)
        txt = txt.replace("disk_io = 'low'", "disk_io = 'medium'")
        open(os.path.join(runs, btag + ".in"), "w", newline="\n").write(txt)
        ok, _ = run_pw(runs, btag + ".in", np=PROD.get("np", 8), npool=1,
                       timeout=3600, keep_outdir=False)
        if not ok:
            print(f"[HATA] {btag} bant kosusu basarisiz", flush=True)
            return res
    kb = parse_bands(boutp)
    if kb:
        # SOC'siz: dolu bant sayisi = 72 elektron / 2 = 36 (cift katman)
        nv = 36
        labels = ["G", "K", "M", "Kd1", "Kd2"]
        for (kv, eigs), lab in zip(kb, labels):
            res[f"vb_{lab}_ev"] = eigs[nv - 1]
            res[f"cb_{lab}_ev"] = eigs[nv]
            if lab in ("G", "K", "M"):
                res[f"vb4_{lab}_ev"] = eigs[nv - 4:nv]
                res[f"cb4_{lab}_ev"] = eigs[nv:nv + 4]
        # K etkin kutleleri (K, Kd1, Kd2 paraboli).
        # DONUSUM (denetimde duzeltildi): hex_cell'in ters orgusunde b1 ve b2
        # arasi 60 derece -> |b1+b2| = 2*(2pi/a). frac (d,d) kaymasi icin
        # |dk_cart| = 2*d*(2pi/a) [QE ciktisiyla sayisal dogrulandi:
        # (1/3,1/3)->(0.3333,0.5774), d=0.0133 -> |dk|=0.0267 (2pi/a)]
        a0 = A0[mat]["a0_ang"]
        import math
        tpa = 2 * math.pi / a0  # 1/A
        dks = [0.0, 0.02 * 2 * tpa, 0.04 * 2 * tpa]
        HBAR2_2M0 = 3.80998212
        try:
            vbp = [res["vb_K_ev"], res["vb_Kd1_ev"], res["vb_Kd2_ev"]]
            cbp = [res["cb_K_ev"], res["cb_Kd1_ev"], res["cb_Kd2_ev"]]
            import numpy as _np
            # K ekstremum: E = E0 + c2 k^2 (lineer terimsiz fit)
            A = _np.column_stack([_np.ones(3), _np.array(dks) ** 2])
            cv = _np.linalg.lstsq(A, _np.array(vbp), rcond=None)[0][1]
            cc = _np.linalg.lstsq(A, _np.array(cbp), rcond=None)[0][1]
            res["mh_K_m0"] = float(-HBAR2_2M0 / cv) if cv else None
            res["me_K_m0"] = float(HBAR2_2M0 / cc) if cc else None
        except KeyError:
            pass
    return res


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    probe = "--probe" in sys.argv
    mats = args if args else ["WSSe", "MoSSe"]
    for mat in mats:
        runs = os.path.join(HERE, f"runs_{mat}")
        os.makedirs(runs, exist_ok=True)
        deltas = DELTAS[:1] if probe else DELTAS
        results = []
        jp = os.path.join(HERE, f"registry_ham_{mat}.json")
        if os.path.exists(jp):
            results = json.load(open(jp))
        # bant adimi da basarili olmali; yoksa nokta yeniden denenir
        done = {r["registry"] for r in results
                if r.get("etot_ry") is not None and "vb_K_ev" in r}
        todo = [(n, d) for n, d in deltas if n not in done]
        print(f"== {mat}: {len(todo)} registry kosusu ==", flush=True)
        for n, d in todo:
            r = one_registry(mat, n, d, runs)
            results = [x for x in results if x["registry"] != n] + [r]
            with open(jp, "w") as f:
                json.dump(results, f, indent=1)
            print(f"{n}: E={r.get('etot_ry')} Ry  dip={r.get('dipole_debye')} D"
                  f"  d={r.get('d_interface_ang')} A  t={r.get('wall_rlx')}",
                  flush=True)
        # CSV ozet
        with open(os.path.join(HERE, f"registry_tarama_{mat}.csv"), "w",
                  newline="") as f:
            w = csv.writer(f)
            w.writerow(["registry", "dx", "dy", "etot_ry", "dipol_debye",
                        "pot_amp_ry", "d_arayuz_A", "vb_K_ev", "cb_K_ev",
                        "vb_G_ev", "cb_G_ev"])
            for r in sorted(results, key=lambda x: x["registry"]):
                w.writerow([r["registry"], r["delta"][0], r["delta"][1],
                            r.get("etot_ry"), r.get("dipole_debye"),
                            r.get("pot_amp_ry"), r.get("d_interface_ang"),
                            r.get("vb_K_ev"), r.get("cb_K_ev"),
                            r.get("vb_G_ev"), r.get("cb_G_ev")])
    print("BITTI", flush=True)


if __name__ == "__main__":
    main()
