"""Run an RDKit reaction for one reagent combination and normalize products.

Input molecules are never modified: RunReactants returns freshly constructed
product molecules, and atom-mapping numbers are removed only on those copies.
"""

from dataclasses import dataclass
from typing import List, Tuple

from rdkit import Chem

from .validation import ParsedReagent

MAX_PRODUCT_TUPLES = 100


@dataclass(frozen=True)
class ComboOutcome:
    reagent_ids: Tuple[str, ...]
    # Deduplicated, normalized valid product-tuple keys (within this combo).
    valid_keys: List[Tuple[str, ...]]
    raw_tuple_count: int
    invalid_tuple_count: int
    truncated: bool
    # Empty string when every produced tuple validates.
    invalid_detail: str

    @property
    def no_match(self) -> bool:
        return self.raw_tuple_count == 0

    @property
    def invalid(self) -> bool:
        return self.raw_tuple_count > 0 and not self.valid_keys


def _strip_mapping(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def normalize_tuple(product_mols) -> Tuple[Tuple[str, ...], str]:
    """Return (sorted canonical isomeric SMILES, error detail).

    Every component is chemically validated (sanitization checks valence,
    aromaticity, charge balance, ...). One failing component fails the tuple.
    """
    smiles_parts: List[str] = []
    errors: List[str] = []
    for component_index, mol in enumerate(product_mols):
        _strip_mapping(mol)
        try:
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        except Exception as exc:  # rdKit raises diverse exception types
            errors.append(f"component {component_index}: {type(exc).__name__}: {exc}")
            continue
        if not smiles:
            errors.append(f"component {component_index}: empty canonical SMILES")
            continue
        smiles_parts.append(smiles)
    if errors:
        return tuple(sorted(smiles_parts)), "; ".join(errors)
    return tuple(sorted(smiles_parts)), ""


def run_combination(reaction, reagents: Tuple[ParsedReagent, ...]) -> ComboOutcome:
    mols = tuple(reagent.mol for reagent in reagents)
    ids = tuple(reagent.id for reagent in reagents)

    tuples = reaction.RunReactants(mols, maxProducts=MAX_PRODUCT_TUPLES + 1)
    raw_count = len(tuples)
    truncated = raw_count > MAX_PRODUCT_TUPLES
    tuples = tuples[:MAX_PRODUCT_TUPLES]

    seen: set = set()
    valid_keys: List[Tuple[str, ...]] = []
    invalid_details: List[str] = []
    invalid_tuple_count = 0
    for product_tuple in tuples:
        key, detail = normalize_tuple(product_tuple)
        if detail:
            invalid_tuple_count += 1
            invalid_details.append(detail)
            continue
        if key not in seen:
            seen.add(key)
            valid_keys.append(key)

    # Keep failure messages stable and bounded.
    detail = "" if valid_keys else ("; ".join(invalid_details[:3]) or "")
    return ComboOutcome(
        reagent_ids=ids,
        valid_keys=valid_keys,
        raw_tuple_count=raw_count,
        invalid_tuple_count=invalid_tuple_count,
        truncated=truncated,
        invalid_detail=detail,
    )
