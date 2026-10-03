# Reaction Product Enumeration API

FastAPI backend that applies a **known reaction template** (reaction SMARTS) to a
reagent library and enumerates every chemically valid product tuple with
provenance (the reagent-ID combination that produced it). It does **not**
predict yields, select the "best" product, or merge stereoisomers.

Built with FastAPI 0.115.12 and RDKit 2025.03.6. Use `.venv/bin/python`.

## Run

```bash
.venv/bin/uvicorn app.main:app --host 127.0.0.1 --port 8000
```

Interactive docs: http://127.0.0.1:8000/docs

## Request

`POST /enumerate`

- `reaction_smarts`: a reaction SMARTS with 1–3 reactant positions and one or
  more product components.
- `reagents`: outer list ordered by reactant position; each inner list holds
  the reagents for that position. Every reagent has a unique `id` within its
  position and a `smiles` string parsed by RDKit.

See [`examples/request_amide.json`](examples/request_amide.json).

```bash
curl -s http://127.0.0.1:8000/enumerate \
  -H 'Content-Type: application/json' \
  --data @examples/request_amide.json | .venv/bin/python -m json.tool
```

## Response

- `total_combinations` / `processed_combinations`: size of the Cartesian
  product (input order, leftmost position varies slowest).
- `no_match_combinations`: combinations where the reaction produced no tuple.
- `products`: each record has:
  - `product_tuple`: the **complete** product tuple as sorted canonical
    isomeric SMILES (order-independent across components, multiplicity kept).
    Stereochemistry, isotopes and formal charges are retained; atom-map numbers
    are stripped. Distinct stereoisomers and distinct tuples are not merged.
  - `sources`: every reagent-ID combination producing that tuple; each source
    is listed once, in first-seen order.
- `failures`: per-combination records distinguishing:
  - `no_match` — the template did not match the combination;
  - `invalid_product` — the template matched but every tuple failed chemical
    validation (e.g. illegal valence); the RDKit error is in `detail`;
  - `invalid_product_tuples` — some tuples failed validation and were dropped
    while valid tuples from the same combination are still reported;
  - `incomplete` — the combination produced more tuples than the limit.
- `complete`: `false` if any combination exceeded the per-combo limit.

Illegal tuples never appear in `products`; they are recorded only in
`failures`. Input reagent molecules are never modified (RDKit returns fresh
product molecules; mapping numbers are removed only from those copies).

## Validation and limits

- The whole batch is rejected with HTTP 400 and positioned error details when
  the reaction SMARTS cannot be parsed, the number of reagent lists does not
  match the template's reactant count, a reagent SMILES is invalid, a reagent
  ID is duplicated within one position, or a position list is empty.
- **500** reagent combinations maximum; requests whose Cartesian product is
  larger are rejected before any reaction runs.
- **100** product tuples maximum per combination; tuples beyond the limit are
  not silently dropped — the combination and overall response are flagged
  `complete: false` with an `incomplete` failure.
- Enumeration, product and source ordering are deterministic; repeating a
  request returns identical ordering.

## Project layout

- `app/validation.py` — request parsing, RDKit validation, limit checks.
- `app/products.py` — RDKit reaction execution and product normalization.
- `app/aggregation.py` — tuple dedup and source aggregation.
- `app/enumeration.py` — orchestration.
- `app/main.py` — HTTP layer.
- `tests/test_enumeration.py` — test suite.

## Tests

```bash
.venv/bin/python -m pytest tests/ -q
.venv/bin/python -m compileall -q app tests
```
