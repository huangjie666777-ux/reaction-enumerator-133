"""Input validation: reaction SMARTS and reagent lists, with located errors.

Any invalid input rejects the whole batch (HTTP 400); nothing is enumerated.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .models import (
    MAX_COMBINATIONS,
    MAX_REACTANT_POSITIONS,
    MIN_REACTANT_POSITIONS,
    ReagentIn,
)


class BatchValidationError(ValueError):
    """One or more input errors; the entire request must be rejected."""

    def __init__(self, details: list[dict[str, str]]):
        self.details = details
        super().__init__("; ".join(d["message"] for d in details))


def _err(location: str, message: str) -> dict[str, str]:
    return {"location": location, "message": message}


def parse_reaction(reaction_smarts: str) -> rdChemReactions.ChemicalReaction:
    """Parse and sanitize the reaction SMARTS; raise BatchValidationError on failure."""
    rxn = None
    try:
        rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    except Exception as exc:  # ChemicalReactionParserException, etc.
        raise BatchValidationError(
            [_err("reaction_smarts", f"Cannot parse reaction SMARTS: {exc}")]
        ) from exc
    if rxn is None:
        raise BatchValidationError(
            [_err("reaction_smarts", "Cannot parse reaction SMARTS.")]
        )
    try:
        rdChemReactions.SanitizeRxn(rxn)
    except Exception as exc:
        raise BatchValidationError(
            [_err("reaction_smarts", f"Invalid reaction SMARTS: {exc}")]
        ) from exc
    n_warnings, n_errors = rxn.Validate(0)
    if n_errors:
        raise BatchValidationError(
            [_err("reaction_smarts", f"Reaction validation reported {n_errors} error(s).")]
        )
    n_reactants = rxn.GetNumReactantTemplates()
    if not (MIN_REACTANT_POSITIONS <= n_reactants <= MAX_REACTANT_POSITIONS):
        raise BatchValidationError(
            [
                _err(
                    "reaction_smarts",
                    f"Reaction must have {MIN_REACTANT_POSITIONS}-{MAX_REACTANT_POSITIONS} "
                    f"reactant templates, found {n_reactants}.",
                )
            ]
        )
    if rxn.GetNumProductTemplates() < 1:
        raise BatchValidationError(
            [_err("reaction_smarts", "Reaction must declare at least one product template.")]
        )
    return rxn


def validate_request(
    reaction_smarts: str,
    reactant_positions: list[list[ReagentIn]],
) -> tuple[rdChemReactions.ChemicalReaction, list[list[tuple[str, Chem.Mol]]]]:
    """Validate the full batch and return (reaction, parsed reagents per position).

    All errors are collected so the client can fix every located problem at once,
    but the batch itself is rejected as a whole.
    """
    errors: list[dict[str, str]] = []
    rxn = None
    try:
        rxn = parse_reaction(reaction_smarts)
    except BatchValidationError as exc:
        errors.extend(exc.details)

    if not reactant_positions:
        errors.append(_err("reactants", "At least one reactant position is required."))
    if len(reactant_positions) > MAX_REACTANT_POSITIONS:
        errors.append(
            _err(
                "reactants",
                f"At most {MAX_REACTANT_POSITIONS} positions are accepted.",
            )
        )

    if rxn is not None and reactant_positions and len(reactant_positions) != rxn.GetNumReactantTemplates():
        errors.append(
            _err(
                "reactants",
                f"Expected {rxn.GetNumReactantTemplates()} reactant position(s) per template, "
                f"received {len(reactant_positions)}.",
            )
        )

    parsed_positions: list[list[tuple[str, Chem.Mol]]] = []
    for pos_idx, reagents in enumerate(reactant_positions):
        if not reagents:
            errors.append(_err(f"reactants[{pos_idx}]", "Reactant position must not be empty."))
            parsed_positions.append([])
            continue
        seen_ids: set[str] = set()
        parsed: list[tuple[str, Chem.Mol]] = []
        for reagent_idx, reagent in enumerate(reagents):
            loc = f"reactants[{pos_idx}][{reagent_idx}]"
            if reagent.id in seen_ids:
                errors.append(_err(f"{loc}.id", f"Duplicate reagent id {reagent.id!r} within position {pos_idx}."))
            else:
                seen_ids.add(reagent.id)
            mol = Chem.MolFromSmiles(reagent.smiles)
            if mol is None:
                errors.append(
                    _err(f"{loc}.smiles", f"Cannot parse SMILES for reagent {reagent.id!r}.")
                )
                continue
            parsed.append((reagent.id, mol))
        parsed_positions.append(parsed)

    if errors:
        raise BatchValidationError(errors)

    total_combinations = 1
    for reagents in parsed_positions:
        total_combinations *= len(reagents)
    if total_combinations > MAX_COMBINATIONS:
        raise BatchValidationError(
            [
                _err(
                    "reactants",
                    f"Cartesian product has {total_combinations} combinations, "
                    f"limit is {MAX_COMBINATIONS}; request rejected before enumeration.",
                )
            ]
        )

    assert rxn is not None
    return rxn, parsed_positions
