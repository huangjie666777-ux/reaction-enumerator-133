"""Aggregate per-combination outcomes into unique products and source lists."""

from collections import OrderedDict
from typing import Dict, List

from .products import ComboOutcome
from .schemas import FailureRecord, ProductRecord


def aggregate(
    outcomes: List[ComboOutcome],
) -> Dict[str, object]:
    """Build the response payload.

    Product keys are sorted product tuples: order-independent across tuple
    components while preserving component multiplicity. Sources never repeat.
    """
    products: "OrderedDict[tuple, List[List[str]]]" = OrderedDict()
    failures: List[FailureRecord] = []
    no_match = 0
    complete = True

    for outcome in outcomes:
        ids = list(outcome.reagent_ids)
        if outcome.truncated:
            complete = False
        if outcome.no_match:
            no_match += 1
            failures.append(
                FailureRecord(
                    reagent_ids=ids,
                    reason="no_match",
                    truncated_product_tuples=outcome.truncated,
                )
            )
            continue
        if outcome.truncated:
            failures.append(
                FailureRecord(
                    reagent_ids=ids,
                    reason="incomplete",
                    detail="more than 100 product tuples; reporting is truncated",
                    truncated_product_tuples=True,
                )
            )
        if outcome.invalid_tuple_count and not outcome.valid_keys:
            failures.append(
                FailureRecord(
                    reagent_ids=ids,
                    reason="invalid_product",
                    detail=outcome.invalid_detail,
                    truncated_product_tuples=outcome.truncated,
                )
            )
        elif outcome.invalid_tuple_count:
            # Valid products exist; invalid tuples are still recorded here and
            # never mixed into the product list.
            failures.append(
                FailureRecord(
                    reagent_ids=ids,
                    reason="invalid_product_tuples",
                    detail=(
                        f"{outcome.invalid_tuple_count} invalid tuple(s) dropped; "
                        f"{outcome.invalid_detail}"
                    ),
                    truncated_product_tuples=outcome.truncated,
                )
            )
        for key in outcome.valid_keys:
            products.setdefault(key, [])
            if ids not in products[key]:
                products[key].append(ids)

    product_records = [
        ProductRecord(product_tuple=list(key), sources=sources)
        for key, sources in products.items()
    ]
    return {
        "total_combinations": len(outcomes),
        "processed_combinations": len(outcomes),
        "no_match_combinations": no_match,
        "complete": complete,
        "products": product_records,
        "failures": failures,
    }
