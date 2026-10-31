"""QupKake smoke test: predict micro-pKa for acetic acid and phenol."""

import os
import sys
import warnings

os.environ["XTBPATH"] = r"D:\GUASSIAN\xtb-6.7.1\bin\xtb.exe"
warnings.filterwarnings("ignore")

import tempfile  # noqa: E402

import torch  # noqa: E402

torch.set_num_threads(2)

from rdkit import Chem  # noqa: E402
from rdkit.Chem import PandasTools  # noqa: E402

from qupkake.predict import run_prediction_pipeline  # noqa: E402


def predict(smiles: str, tag: str):
    tmp = tempfile.mkdtemp(prefix=f"qupkake_{tag}_")
    os.makedirs(os.path.join(tmp, "raw"), exist_ok=True)
    os.makedirs(os.path.join(tmp, "output"), exist_ok=True)
    mol = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(Chem.Mol(mol))
    from rdkit.Chem import AllChem

    AllChem.EmbedMolecule(m, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(m)
    m = Chem.RemoveHs(m)
    m.SetProp("_Name", tag)
    w = Chem.SDWriter(os.path.join(tmp, "raw", f"{tag}.sdf"))
    w.write(m)
    w.close()
    run_prediction_pipeline(
        root=tmp, filename=f"{tag}.sdf", tautomerize=False,
        name_col="_Name", mol_col="ROMol", mp=False, output=f"{tag}_out.sdf",
    )
    out = os.path.join(tmp, "output", f"{tag}_out.sdf")
    rows = []
    for m2 in Chem.SDMolSupplier(out, removeHs=False):
        if m2 is None:
            continue
        rows.append((m2.GetProp("idx"), m2.GetProp("pka_type"),
                     float(m2.GetProp("pka"))))
    return rows


for name, smi, expect in [("acetic acid", "CC(=O)O", "~4.7"),
                          ("phenol", "c1ccccc1O", "~10.0"),
                          ("imidazole", "c1ncc[nH]1", "~7.0 (conj acid)")]:
    try:
        rows = predict(smi, name.replace(" ", "_"))
        print(f"{name:12s} (expect {expect}):", rows, flush=True)
    except Exception as e:
        import traceback

        traceback.print_exc()
        print(f"{name}: FAILED {type(e).__name__}: {e}", flush=True)
