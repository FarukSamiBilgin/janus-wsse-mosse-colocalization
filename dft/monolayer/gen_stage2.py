# -*- coding: utf-8 -*-
"""Asama 1b uretici: vc-relax ciktisindan gevsemis geometriyi al,
SOC (full-rel) scf + bant + Rashba-ince-tarama girdilerini yaz.
Ayrica skaler ciktidan dipol/orgu sabitini parametre sayfasina isle."""
import re, json, os, sys

D = os.path.dirname(os.path.abspath(__file__))
BOHR = 0.529177210903

MAT = {
    'mosse': dict(M='Mo', MPS='Mo.rel-pbe-spn-kjpaw_psl.1.0.0.UPF', Mm=95.95),
    'wsse':  dict(M='W',  MPS='W.rel-pbe-spn-kjpaw_psl.1.0.0.UPF',  Mm=183.84),
}

def parse_vcrelax(path):
    txt = open(path, encoding='utf-8', errors='ignore').read()
    # son CELL_PARAMETERS blogu (alat cinsinden) + alat
    alats = re.findall(r'lattice parameter \(alat\)\s*=\s*([\d.]+)\s*a\.u\.', txt)
    alat0 = float(alats[0]) * BOHR   # baslangic alat (Angstrom)
    cells = re.findall(r'CELL_PARAMETERS \(alat=\s*([\d.]+)\)\s*\n'
                       r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n'
                       r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n'
                       r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', txt)
    if cells:
        c = cells[-1]
        alat = float(c[0]) * BOHR
        a_vec = [float(c[1])*alat, float(c[2])*alat, float(c[3])*alat]
        a_len = (a_vec[0]**2 + a_vec[1]**2)**0.5
        c_len = float(c[9-3+2])*alat  # a3_z: index 7,8,9 -> c[7],c[8],c[9]... duzelt
        c_len = float(c[9])*alat
    else:
        a_len = alat0; c_len = 20.0
    # son ATOMIC_POSITIONS (crystal)
    pos_blocks = re.findall(r'ATOMIC_POSITIONS \(crystal\)\s*\n((?:\s*\w+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)', txt)
    atoms = []
    if pos_blocks:
        for line in pos_blocks[-1].strip().splitlines():
            p = line.split()
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    # dipol (son deger, e*Angstrom'a cevrilebilir; QE 'Debye' ve a.u. yazar)
    dips = re.findall(r'dipole moment.*?=\s*([-\d.]+)\s*Debye', txt)
    dip_debye = float(dips[-1]) if dips else None
    # toplam enerji
    enes = re.findall(r'!\s*total energy\s*=\s*([-\d.]+)\s*Ry', txt)
    return dict(a=a_len, c=c_len, atoms=atoms, dip_debye=dip_debye,
                E_Ry=float(enes[-1]) if enes else None)

HEADER = """&control
  calculation = '{calc}'
  prefix = '{pref}'
  outdir = './tmp_{pref}'
  pseudo_dir = './pseudo'
  tefield = .true.
  dipfield = .true.
  {extra_c}
/
&system
  ibrav = 4
  a = {a:.6f}
  c = {c:.6f}
  nat = 3
  ntyp = 3
  ecutwfc = 60.0
  ecutrho = 480.0
  occupations = 'fixed'
  edir = 3
  emaxpos = 0.05
  eopreg = 0.03
  eamp = 0.0
  noncolin = .true.
  lspinorb = .true.
  {extra_s}
/
&electrons
  conv_thr = {conv}
  mixing_beta = 0.3
  {extra_e}
/
ATOMIC_SPECIES
{M:3s} {Mm:.2f}  {MPS}
S   32.06  S.rel-pbe-n-kjpaw_psl.1.0.0.UPF
Se  78.97  Se.rel-pbe-dn-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS (crystal)
{positions}
"""

def kpath_bands(a_len):
    """Gamma-K-M-Gamma + Gamma cevresi ince cizgiler (Rashba fiti icin)."""
    lines = []
    # kaba yol (crystal_b)
    pts = [('0.0 0.0 0.0', 30), ('0.3333333 0.3333333 0.0', 15),
           ('0.5 0.0 0.0', 20), ('0.0 0.0 0.0', 1)]
    return pts

def write_inputs(m, g):
    info = MAT[m]
    pos = '\n'.join('%-3s %12.9f %12.9f %12.9f' % a for a in g['atoms'])
    pref = m + '_soc'
    # SOC scf
    scf = HEADER.format(calc='scf', pref=pref, a=g['a'], c=g['c'],
                        conv='1.0d-9', extra_c='', extra_s='', extra_e='',
                        M=info['M'], Mm=info['Mm'], MPS=info['MPS'],
                        positions=pos)
    scf += "K_POINTS (automatic)\n12 12 1 0 0 0\n"
    open(os.path.join(D, f'{m}_soc_scf.in'), 'w').write(scf)
    # SOC bantlar: kaba yol
    b = HEADER.format(calc='bands', pref=pref, a=g['a'], c=g['c'],
                      conv='1.0d-9', extra_c='', extra_s='nbnd = 40',
                      extra_e='diago_full_acc = .true.',
                      M=info['M'], Mm=info['Mm'], MPS=info['MPS'], positions=pos)
    b += ("K_POINTS (crystal_b)\n4\n"
          "0.000000 0.000000 0.0 30\n"
          "0.333333 0.333333 0.0 15\n"
          "0.500000 0.000000 0.0 20\n"
          "0.000000 0.000000 0.0 1\n")
    open(os.path.join(D, f'{m}_soc_bands.in'), 'w').write(b)
    # Rashba ince taramasi: Gamma cevresi +-0.08 (2pi/a birimi) G-K ve G-M yonleri
    ks = []
    n = 17
    for i in range(n):
        t = -0.08 + 0.16*i/(n-1)
        ks.append((t/3**0.5, t/3, 0.0))     # kabaca G-K dogrultusuna yakin
    for i in range(n):
        t = -0.08 + 0.16*i/(n-1)
        ks.append((t, 0.0, 0.0))            # G-M (kartezyen 2pi/a)
    r = HEADER.format(calc='bands', pref=pref, a=g['a'], c=g['c'],
                      conv='1.0d-9', extra_s='nbnd = 40', extra_c='',
                      extra_e='diago_full_acc = .true.',
                      M=info['M'], Mm=info['Mm'], MPS=info['MPS'], positions=pos)
    r += 'K_POINTS (tpiba_b)\n%d\n' % len(ks)
    for k in ks:
        r += '%.6f %.6f %.6f 1\n' % k
    open(os.path.join(D, f'{m}_soc_rashba.in'), 'w').write(r)

if __name__ == '__main__':
    sonuc = {}
    for m in ('mosse', 'wsse'):
        g = parse_vcrelax(os.path.join(D, f'{m}_vcrelax.out'))
        sonuc[m] = dict(a_Angstrom=g['a'], dip_Debye=g['dip_debye'], E_Ry=g['E_Ry'],
                        atoms=g['atoms'])
        write_inputs(m, g)
        print(m, 'a = %.4f A' % g['a'], 'dipol = %s D' % g['dip_debye'])
    json.dump(sonuc, open(os.path.join(D, 'monolayer_geo.json'), 'w'), indent=2)
    print('SOC girdileri yazildi.')
