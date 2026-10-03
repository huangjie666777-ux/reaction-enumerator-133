# Reaction Product Enumeration Backend

Applies a **known reaction SMARTS template** to a library of reagents and
enumerates every product it can form, with full source provenance. It does
**not** predict yields, rates, or selectivity; each candidate is simply a
structure the template can produce.

Built with FastAPI 0.115.12 and RDKit 2025.03.6. Run everything with
`.venv/bin/python`.

## Layout

- `app/models.py` — request/response schemas and limits
- `app/validation.py` — reaction SMARTS parsing, reagent SMILES parsing, located errors
- `app/reaction_runner.py` — cartesian-product enumeration and RDKit reaction execution
- `app/normalization.py` — full sanitization, map removal, canonical isomeric SMILES
- `app/aggregation.py` — unique-product / source aggregation
- `app/main.py` — FastAPI app, `POST /enumerate`, `GET /health`
- `tests/test_enumeration.py` — 18 self-contained pytest cases

## Request

`POST /enumerate`

```json
{
  "reaction_smarts": "[N:1].[C:2](=[O:3])Cl>>[N:1][C:2](=[O:3])",
  "reactants": [
    [{"id": "methylamine", "smiles": "CN"}],
    [{"id": "acetyl-chloride", "smiles": "CC(=O)Cl"}]
  ]
}
```

- `reactants` is grouped **by reactant position**, in the same order as the
  dot-separated templates on the left side of the SMARTS.
- Reagent `id` must be unique within its position; `smiles` is parsed by RDKit.
- Reactions have 1-3 reactant templates and at least one product template.
- The number of positions must equal the template reactant count.

Invalid input rejects the **entire batch** with HTTP 400 and an error body
whose `details` entries carry `location` (either `reaction_smarts` or e.g.
`reactants[0][2].smiles`) and `message`.

## Semantics

- Reagent combinations are enumerated in input order (stable cartesian product).
- RDKit runs the template at **all** sites; every product tuple is kept.
- A combination with no reacting site counts as `no_match_combinations`.
- Every product molecule is fully sanitized. Illegal valence or other chemical
  invalidity is recorded per combination under `failures`
  (`type: invalid_product`) and never mixed into valid results.
- Products are canonical **isomeric** SMILES after stripping atom-map numbers,
  preserving stereochemistry, isotopes and formal charges.
- The unique key is the **sorted full product tuple**: component order is
  ignored while component multiplicity is retained; stereoisomers are never
  merged; multi-component products keep every component (no largest-fragment rule).
- Each unique product lists every source reagent-id combination once.
- Input reagent molecules are copied before reaction and never modified.

## Response

```json
{
  "total_combinations": 3,
  "processed_combinations": 3,
  "no_match_combinations": 1,
  "complete": true,
  "products": [
    {"smiles": ["CNC(C)=O"], "sources": [["methylamine", "acetyl-chloride"]], "source_count": 1}
  ],
  "failures": [],
  "limits": {"max_combinations": 500, "max_product_tuples_per_combination": 100}
}
```

## Limits

- At most **500** reagent combinations per request (cartesian product size).
  Requests exceeding this are rejected with HTTP 400 **before** any work runs.
- At most **100 product tuples per combination**. If a combination produces
  more, the first 100 are kept, a `too_many_product_tuples` failure is
  reported for that combination, and the top-level `complete` flag is `false`;
  results are never silently truncated.
- Repeated identical requests return identical products and source ordering.

## Run

```bash
.venv/bin/python -m uvicorn app.main:app --host 127.0.0.1 --port 8000
curl -s http://127.0.0.1:8000/health
curl -s http://127.0.0.1:8000/enumerate -H 'Content-Type: application/json' \
  --data @examples/payload.json
```

## Test

```bash
.venv/bin/python -m pytest tests/ -q
.venv/bin/python -m compileall app tests
```
