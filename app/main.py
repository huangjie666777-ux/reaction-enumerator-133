"""FastAPI application for reaction-template product enumeration."""
from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from .aggregation import aggregate
from .models import EnumerationRequest
from .validation import BatchValidationError, validate_request

app = FastAPI(
    title="Reaction Product Enumerator",
    version="1.0.0",
    description="Apply known reaction SMARTS templates to reagent libraries. "
    "Generates sourced candidate structures; does not predict yields.",
)


@app.exception_handler(BatchValidationError)
async def batch_validation_handler(_request: Request, exc: BatchValidationError) -> JSONResponse:
    return JSONResponse(status_code=400, content={"error": "invalid_input", "details": exc.details})


@app.exception_handler(RequestValidationError)
async def request_validation_handler(_request: Request, exc: RequestValidationError) -> JSONResponse:
    details = [
        {
            "location": ".".join(str(part) for part in error["loc"] if part != "body"),
            "message": error["msg"],
        }
        for error in exc.errors()
    ]
    return JSONResponse(status_code=400, content={"error": "invalid_input", "details": details})


@app.get("/health")
async def health() -> dict:
    return {"status": "ok"}


@app.post("/enumerate")
async def enumerate_products(request: EnumerationRequest) -> dict:
    rxn, parsed_positions = validate_request(request.reaction_smarts, request.reactants)
    return aggregate(rxn, parsed_positions)
