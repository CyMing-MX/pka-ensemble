# pka_replica v1.0

Open-source aqueous pKa prediction for organic molecules: site-level (micro) and
molecule-level (macro) routes, microstate ensemble thermodynamics, a local web
UI, batch processing and xlsx reporting.
No Python required — unzip and double-click `pka-server.exe`.

## Downloads — which one?

| File | Size | Contains |
|---|---|---|
| `pka_replica_v1.0_lite.zip` | ~199 MB | sites (SMARTS fragment) + xtb (GFN2-xTB/ALPB) routes |
| `pka_replica_v1.0_full.zip` | ~473 MB | lite + `_ai/` — adds the AIMNet2 pre-optimization stage (wB97M-D3 neural potential), identical to the developer setup |
| `pka_replica_v1.0_ai-addon.zip` | ~274 MB | just the `_ai/` folder — extract into an existing lite install to upgrade in place |

Verify integrity with `SHA256SUMS.txt`.

## Method

Microstate ensemble thermodynamics following the Uni-pKa (JACS Au) /
Starling (Wagen 2025) formalism: SMARTS-template microstate enumeration →
per-microstate free energies (sites: fragment constants + Coulomb coupling;
xtb: geometry reuse + optional AIMNet2 pre-optimization → GFN2-xTB gas-phase
optimization → ALPB-water single point, with per-functional-group linear
calibration) → Boltzmann partitioning → macro pKa / pI / logD curves.

## Accuracy (experimental references only — no training data)

- pKahub benchmark (Sipos-Szabó et al., JCIM 2026, monoprotic):
  best-match MAE **1.08** (xtb route, n=46 sample) / 1.72 (sites, n=6,084);
  reference: MolGpka 0.75, QupKake 0.79, pKaSolver 0.89, Chemaxon 0.94
- 32 oncology drugs vs Chemicalize: pKa-B MAE **1.17**, pKa-A 1.72, pI 1.45
- 11 perfluoro acids (PFAS) vs Chemicalize: MAE **0.63** (all perfluoro
  carboxylic acids within 0.5)
- single-site spot checks: acetic 4.42 (4.76), pyridine 5.23 (5.23),
  imidazole 6.93 (6.95), ibuprofen 4.36 (4.43)

## Quick start

1. Extract the zip to a local disk path (avoid non-ASCII/spaces)
2. Double-click `pka-server.exe` — a console opens and the browser pops up at `http://127.0.0.1:8765`
3. Type a SMILES / CAS / compound name and hit predict

CLI: `pka-cli.exe "58-08-2"`, `pka-cli.exe "caffeine" --micro`, `pka-cli.exe "CC(=O)O" --json out.json`

See `操作手册.md` (user manual, Chinese) inside each package for the full guide.

## Notes

- PubChem-backed features need internet; everything else runs fully offline
- The AIMNet2 stage (full / addon) compiles warp kernels on first run (~10–30 s)
- Known limits: zwitterions/polybasic coupling is approximate; no RRHO
  corrections; single-conformer ensemble; empirical calibration is only
  trustworthy inside its reference chemistry space
- Third-party licenses (xtb LGPL-3.0, RDKit BSD, torch BSD, …) are in `LICENSES/`
  inside each package; this project itself is MIT
