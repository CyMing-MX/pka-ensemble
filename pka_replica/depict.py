"""2D structure depiction helpers (RDKit SVG), shared by the web server and
the report/export modules."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

HL_COLOR = (0.99, 0.61, 0.16)  # warm orange highlight


def depict_mol(mol: Chem.Mol, width: int = 240, height: int = 170,
               highlight: list[int] | None = None) -> str | None:
    """2D depiction as SVG with a transparent background (default light-atom
    palette suits the light UI)."""
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    m = Chem.Mol(mol)
    try:
        AllChem.Compute2DCoords(m)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False  # transparent -> blends with the card
        opts.padding = 0.06
        hl = [i for i in (highlight or []) if 0 <= i < m.GetNumAtoms()]
        if hl:
            colors = {i: HL_COLOR for i in hl}
            drawer.DrawMolecule(m, highlightAtoms=hl, highlightAtomColors=colors)
        else:
            drawer.DrawMolecule(m)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def depict_smiles(smiles: str, **kw) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return depict_mol(mol, **kw)
