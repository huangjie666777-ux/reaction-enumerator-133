"""Input validation and RDKit parsing. All invalid input rejects the whole batch."""

from dataclasses import dataclass
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .schemas import EnumerationRequest, ValidationErrorItem

MIN_REACTANTS = 1
MAX_REACTANTS = 3
MAX_COMBINATIONS = 500


class RequestValidationError(Exception):
    def __init__(self, details: List[ValidationErrorItem]):
        self.details = details
        super().__init__("; ".join(d.message for d in details))


@dataclass(frozen=True)
class ParsedReagent:
    id: str
    smiles: str  # original input SMILES
    mol: Chem.Mol  # parse of the input; never mutated downstream


@dataclass(frozen=True)
class ValidatedRequest:
    reaction: rdChemReactions.ChemicalReaction
    slots: List[List[ParsedReagent]]


def _err(location, message, position=-1, index=-1, reagent_id=""):
    return ValidationErrorItem(
        location=location,
        position=position,
        index=index,
        reagent_id=reagent_id,
        message=message,
    )


def validate_request(payload: EnumerationRequest) -> ValidatedRequest:
    errors: List[ValidationErrorItem] = []

    try:
        reaction = rdChemReactions.ReactionFromSmarts(payload.reaction_smarts)
    except Exception:
        reaction = None
    if reaction is None or reaction.GetNumReactantTemplates() == 0:
        errors.append(
            _err(
                "reaction_smarts",
                "cannot parse reaction SMARTS or it has no reactant side",
            )
        )
        raise RequestValidationError(errors)
    if not (MIN_REACTANTS <= reaction.GetNumReactantTemplates() <= MAX_REACTANTS):
        errors.append(
            _err(
                "reaction_smarts",
                f"reaction declares {reaction.GetNumReactantTemplates()} reactants; "
                f"only {MIN_REACTANTS}..{MAX_REACTANTS} are supported",
            )
        )
        raise RequestValidationError(errors)
    if reaction.GetNumProductTemplates() < 1:
        errors.append(_err("reaction_smarts", "reaction declares no products"))
        raise RequestValidationError(errors)

    n_positions = reaction.GetNumReactantTemplates()
    slots = payload.reagents
    if len(slots) != n_positions:
        errors.append(
            _err(
                "reaction_smarts",
                f"reactant position count mismatch: template needs {n_positions}, "
                f"request supplies {len(slots)} reagent lists",
            )
        )
        raise RequestValidationError(errors)

    parsed_slots: List[List[ParsedReagent]] = []
    total_combinations = 1
    for position, slot in enumerate(slots):
        if not isinstance(slot, list) or len(slot) == 0:
            errors.append(
                _err("reagent", "empty reagent list for reactant position", position)
            )
            continue
        seen_ids = set()
        parsed_slot: List[ParsedReagent] = []
        for index, reagent in enumerate(slot):
            if reagent.id in seen_ids:
                errors.append(
                    _err(
                        "reagent",
                        "duplicate reagent id within one reactant position",
                        position,
                        index,
                        reagent.id,
                    )
                )
                continue
            seen_ids.add(reagent.id)
            mol = Chem.MolFromSmiles(reagent.smiles)
            if mol is None:
                errors.append(
                    _err(
                        "reagent",
                        "cannot parse reagent SMILES with RDKit",
                        position,
                        index,
                        reagent.id,
                    )
                )
                continue
            parsed_slot.append(ParsedReagent(reagent.id, reagent.smiles, mol))
        if parsed_slot:
            total_combinations *= len(parsed_slot)
        parsed_slots.append(parsed_slot)

    if errors:
        raise RequestValidationError(errors)

    if total_combinations > MAX_COMBINATIONS:
        errors.append(
            _err(
                "request",
                f"cartesian product has {total_combinations} combinations; "
                f"limit is {MAX_COMBINATIONS}",
            )
        )
        raise RequestValidationError(errors)

    return ValidatedRequest(reaction=reaction, slots=parsed_slots)


def iter_combinations(
    slots: List[List[ParsedReagent]],
) -> List[Tuple[ParsedReagent, ...]]:
    """Cartesian product in input order (leftmost position varies slowest)."""
    combos: List[Tuple[ParsedReagent, ...]] = [()]
    for slot in slots:
        combos = [prefix + (reagent,) for prefix in combos for reagent in slot]
    return combos
