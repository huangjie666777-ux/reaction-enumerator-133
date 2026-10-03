"""Apply a reaction template to reagent combinations in input order."""
from __future__ import annotations

from itertools import product
from typing import Iterator

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .models import MAX_PRODUCT_TUPLES_PER_COMBINATION
from .normalization import InvalidProductError, canonicalize_product, product_tuple_key


def enumerate_combinations(
    positions: list[list[tuple[str, Chem.Mol]]],
) -> Iterator[tuple[list[str], list[Chem.Mol]]]:
    """Yield (ordered reagent ids, fresh reactant molecule copies) per cartesian product.

    Order follows the input lists (leftmost position slowest in itertools terms),
    giving stable, repeatable enumeration.
    """
    for combo in product(*positions):
        ids = [reagent_id for reagent_id, _ in combo]
        reactants = [Chem.Mol(mol) for _, mol in combo]
        yield ids, reactants


def run_combination(
    rxn: rdChemReactions.ChemicalReaction,
    reactants: list[Chem.Mol],
) -> tuple[str, list[tuple[str, ...]], list[dict]]:
    """Run one reagent combination.

    Returns (status, valid_product_keys, failures) where status is one of
    ``"ok"``, ``"no_match"`` or ``"invalid_product"``.  Failures include the
    product-tuple index and component message; valid tuples from the same
    combination are still returned.
    """
    try:
        raw_tuples = rxn.RunReactants(
            tuple(reactants),
            maxProducts=MAX_PRODUCT_TUPLES_PER_COMBINATION + 1,
        )
    except Exception as exc:
        return (
            "invalid_product",
            [],
            [{"type": "reaction_error", "message": str(exc).splitlines()[0], "product_tuple_index": None}],
        )

    if not raw_tuples:
        return "no_match", [], []

    truncated = len(raw_tuples) > MAX_PRODUCT_TUPLES_PER_COMBINATION
    raw_tuples = raw_tuples[:MAX_PRODUCT_TUPLES_PER_COMBINATION]

    valid_keys: list[tuple[str, ...]] = []
    failures: list[dict] = []
    for tuple_index, product_tuple in enumerate(raw_tuples):
        component_smiles: list[str] = []
        tuple_failed = False
        for component_index, mol in enumerate(product_tuple):
            try:
                component_smiles.append(canonicalize_product(mol))
            except InvalidProductError as exc:
                tuple_failed = True
                failures.append(
                    {
                        "type": "invalid_product",
                        "message": f"Product tuple {tuple_index} component {component_index}: {exc}",
                        "product_tuple_index": tuple_index,
                    }
                )
        if not tuple_failed:
            valid_keys.append(product_tuple_key(component_smiles))

    if truncated:
        failures.append(
            {
                "type": "too_many_product_tuples",
                "message": (
                    f"More than {MAX_PRODUCT_TUPLES_PER_COMBINATION} product tuples were produced; "
                    "only the first ones were kept and this combination is incomplete."
                ),
                "product_tuple_index": None,
            }
        )

    status = "invalid_product" if failures else "ok"
    return status, valid_keys, failures
