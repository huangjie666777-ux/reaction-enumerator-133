"""Product molecule validation and canonicalization.

Products are fully sanitized (invalid valence etc. is detected), atom-map
numbers are stripped, and canonical isomeric SMILES preserve stereochemistry,
isotopes and formal charges.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem


class InvalidProductError(ValueError):
    """A product tuple component failed full chemical validation."""


def canonicalize_product(product: Chem.Mol) -> str:
    """Sanitize, strip mapping numbers and return canonical isomeric SMILES.

    Raises InvalidProductError on illegal valence or other sanitization failure.
    """
    mol = Chem.Mol(product)
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise InvalidProductError(str(exc).splitlines()[0]) from exc
    try:
        AllChem.RemoveAtomMapping(mol)
    except Exception:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    if not smiles:
        raise InvalidProductError("Empty structure after canonicalization.")
    return smiles


def product_tuple_key(component_smiles: list[str]) -> tuple[str, ...]:
    """Order-independent key that keeps component multiplicity.

    Sorting only reorders components; duplicates remain duplicated, and
    distinct stereoisomers have distinct isomeric SMILES and are not merged.
    """
    return tuple(sorted(component_smiles))
