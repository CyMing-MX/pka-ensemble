"""CLI for pka_replica.

Examples:
    pka-cli "CC(=O)O"                       # macro pKa (auto backend)
    pka-cli "58-08-2"                       # CAS number -> PubChem -> pKa
    pka-cli "caffeine" --micro              # compound name works too
    pka-cli "NCC(=O)O" --backend sites      # SMARTS fragment route
    pka-cli "CC(=O)O" --json out.json       # full result JSON
"""

from __future__ import annotations

import argparse
import json
import sys

sys.path.insert(0, ".")

from pka_replica.workflow import run_macropka, run_micropka  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description="pka_replica: aqueous micro & macro pKa prediction (microstate ensemble + GFN2-xTB/ALPB)")
    ap.add_argument("smiles")
    ap.add_argument("--backend", default="auto", choices=["auto", "xtb", "sites"])
    ap.add_argument("--micro", action="store_true", help="site-level micro pKa only")
    ap.add_argument("--min-charge", default="auto")
    ap.add_argument("--max-charge", default="auto")
    ap.add_argument("--min-ph", type=float, default=0.0)
    ap.add_argument("--max-ph", type=float, default=14.0)
    ap.add_argument("--n-confs", type=int, default=2)
    ap.add_argument("--json", dest="json_out", default=None, help="write full result JSON here")
    args = ap.parse_args()

    from rdkit import Chem

    _pubchem_info = None
    query = args.smiles
    if Chem.MolFromSmiles(query) is None:
        from . import pubchem

        try:
            info = pubchem.resolve(query)
        except pubchem.PubChemError as e:
            sys.exit(f"输入既不是有效 SMILES,PubChem 也未能解析: {e}")
        print(f"[PubChem] {query} -> {info['smiles']} (CID {info['cid']})")
        args.smiles = info["smiles"]
        _pubchem_info = info

    macro = None
    if not args.micro:
        macro = run_macropka(
            args.smiles,
            min_pH=args.min_ph,
            max_pH=args.max_ph,
            min_charge=args.min_charge,
            max_charge=args.max_charge,
            backend=args.backend,
            n_confs=args.n_confs,
        )
        if _pubchem_info:
            # keep the identity block in --json output on par with /api/predict
            macro["resolved_info"] = _pubchem_info
        else:
            try:
                from . import pubchem as _pc
                macro["resolved_info"] = _pc.lookup_smiles(args.smiles)
            except Exception:
                pass
    micro = run_micropka(args.smiles, backend=args.backend,
                         site_pka=macro.get("site_pka_full") if macro else None)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"macro": macro, "micro": micro}, fh, indent=1)
        print(f"[written to {args.json_out}]")

    if macro is not None:
        print(f"\ninput          : {macro['object_data']['initial_smiles']}")
        print(f"backend        : {macro['object_data']['energy_backend']}")
        print(f"microstates    : {len(macro['microstates'])}")
        print("\nmacro-pKa values:")
        for p in macro["pKa_values"]:
            print(
                f"  z {p['initial_charge']:+d} -> {p['final_charge']:+d} :"
                f"  pKa = {p['pKa']:.2f}"
            )
        if macro["isoelectric_point"] is not None:
            print(f"isoelectric point pI = {macro['isoelectric_point']}")
        print(f"logP (输入态, Crippen) = {macro['logP_neutral']}")

        props = macro.get("mol_properties", {})
        if props:
            print("physchem       : MW {MW} | TPSA {TPSA} A^2 | HBD {HBD} | HBA {HBA} "
                  "| rotB {rotatable_bonds} | arom {aromatic_rings} | q {formal_charge}".format(**props))

    print("\nmicro-pKa sites:")
    for b in micro["conjugate_bases"]:
        print(f"  acid  atom {b['atom_index']:3d} [{b['template']}]  pKa = {b['pKa']:.2f}")
    for a in micro["conjugate_acids"]:
        print(f"  base  atom {a['atom_index']:3d} [{a['template']}]  pKa = {a['pKa']:.2f}")


if __name__ == "__main__":
    main()
