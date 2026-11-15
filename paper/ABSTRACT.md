# Extended abstract

**Status: Manuscript in preparation — in revision.** A centre-of-mass mass error was found in an internal audit; quantities depending on it (addressing separations, zero-point energies, emitter band widths) are being recomputed and are not reported here.

The phases of the moiré potential decide whether electrons and holes in a twisted bilayer collect on the same superlattice site (a bright, self-addressed emitter array) or on opposite sites (a charge-separating texture). For polar, R-stacked Janus bilayers this question has rested on a single published parameter table [Angeli et al., Phys. Rev. B 106, 235159 (2022); Erratum 109, 199902 (2024)].

For twisted WSSe and MoSSe homobilayers in the polar (parallel-dipole) stacking I compute, with dipole-corrected PBE+D3 density-functional theory (Quantum ESPRESSO), the local band edges at 12 interlayer registries per material, align them across registries with a semicore-eigenvalue method, and fit a two-shell C₃-closed Fourier model. The fitted K-valley first-star amplitudes (4.0–5.2 meV) are close to the published table. The K–K valence/conduction phase separation is 174° (WSSe; table 172°) and 178° (MoSSe). The tabulated MoSSe valence row describes the Γ-valley hole, so the MoSSe comparison is made on the site assignment rather than on phase values.

A convention-free reading of the high-symmetry registries shows that both carriers localize on the same B-type site in both materials, whichever valence valley (Γ or K) hosts the hole; this confirms the table's localization-site column and excludes the opposite-site assignment implied by its tabulated phases. The local K gap on that site is 103.8 meV (WSSe) and 100.2 meV (MoSSe) smaller than on the opposite site, a reference-free number. With spin–orbit coupling, WSSe has its valence maximum at K, 145 meV above Γ; the corresponding MoSSe check is in progress.

Two reusable methodological safeguards are documented: a C₃-closed first-star triplet in an ibrav-4 basis (verified by a ring-closure test), and slab eigenvalue alignment via semicore states.
