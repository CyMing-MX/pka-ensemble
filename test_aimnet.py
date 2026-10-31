"""AIMNet2 smoke test v2: ASE-based geometry optimization."""
import sys
import time

sys.path.insert(0, r"D:\zcode_rea\pka_replica")
import torch

torch.set_num_threads(4)
from aimnet.calculators import AIMNet2ASE  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.optimize import BFGS  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402


def main():
    t0 = time.time()
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)[O-]"))
    AllChem.ETKDGv3().randomSeed = 42
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    atoms = Atoms(numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                  positions=conf.GetPositions())
    atoms.calc = AIMNet2ASE("aimnet2", charge=-1)
    print(f"model ready in {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    e0 = atoms.get_potential_energy()
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=0.05, steps=300)
    e1 = atoms.get_potential_energy()
    print(f"opt {time.time()-t0:.1f}s: E0={e0:.4f} eV -> E={e1:.4f} eV", flush=True)
    print("AIMNet2 smoke test OK", flush=True)


if __name__ == "__main__":
    main()
