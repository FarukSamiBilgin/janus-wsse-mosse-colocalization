#!/bin/bash
# Kol A DFT: ortam + pseudo indirme (skaler + full-relativistic PAW, pslibrary)
set -e
mkdir -p ~/kolA_dft/pseudo
cd ~/kolA_dft/pseudo
# kolB'den mevcutlari kopyala
cp -n ~/kolB_dft/pseudo/*.UPF . 2>/dev/null || true
BASE="https://pseudopotentials.quantum-espresso.org/upf_files"
for f in \
  Se.pbe-dn-kjpaw_psl.1.0.0.UPF \
  W.pbe-spn-kjpaw_psl.1.0.0.UPF \
  Mo.rel-pbe-spn-kjpaw_psl.1.0.0.UPF \
  S.rel-pbe-n-kjpaw_psl.1.0.0.UPF \
  Se.rel-pbe-dn-kjpaw_psl.1.0.0.UPF \
  W.rel-pbe-spn-kjpaw_psl.1.0.0.UPF ; do
  if [ ! -f "$f" ]; then
    echo "indiriliyor: $f"
    curl -sfLo "$f" "$BASE/$f" || echo "HATA: $f indirilemedi"
  fi
done
ls -la
