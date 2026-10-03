"""End-to-end enumeration orchestration."""

from typing import Dict

from .aggregation import aggregate
from .products import run_combination
from .schemas import EnumerationRequest
from .validation import iter_combinations, validate_request


def enumerate_request(payload: EnumerationRequest) -> Dict[str, object]:
    validated = validate_request(payload)
    combinations = iter_combinations(validated.slots)
    outcomes = [run_combination(validated.reaction, combo) for combo in combinations]
    return aggregate(outcomes)
