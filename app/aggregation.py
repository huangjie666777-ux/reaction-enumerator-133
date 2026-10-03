"""Orchestrate enumeration and aggregate unique products with sources."""
from __future__ import annotations

from collections import OrderedDict

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .models import MAX_COMBINATIONS, MAX_PRODUCT_TUPLES_PER_COMBINATION
from .reaction_runner import enumerate_combinations, run_combination


def aggregate(
    rxn: rdChemReactions.ChemicalReaction,
    positions: list[list[tuple[str, Chem.Mol]]],
) -> dict:
    """Enumerate all combinations and build the unique-product/source report."""
    unique_products: "OrderedDict[tuple[str, ...], list[list[str]]]" = OrderedDict()
    failures: list[dict] = []
    total = 0
    no_match = 0
    incomplete = False

    for ids, reactants in enumerate_combinations(positions):
        total += 1
        status, keys, combo_failures = run_combination(rxn, reactants)
        if status == "no_match":
            no_match += 1
            continue
        for key in keys:
            if key not in unique_products:
                unique_products[key] = []
            sources = unique_products[key]
            if ids not in sources:
                sources.append(list(ids))
        for failure in combo_failures:
            failures.append({"reactant_ids": list(ids), **failure})
            if failure["type"] == "too_many_product_tuples":
                incomplete = True

    products = [
        {
            "smiles": list(key),
            "sources": sources,
            "source_count": len(sources),
        }
        for key, sources in unique_products.items()
    ]

    return {
        "total_combinations": total,
        "processed_combinations": total,
        "no_match_combinations": no_match,
        "complete": not incomplete,
        "products": products,
        "failures": failures,
        "limits": {
            "max_combinations": MAX_COMBINATIONS,
            "max_product_tuples_per_combination": MAX_PRODUCT_TUPLES_PER_COMBINATION,
        },
    }
