"""FastAPI application exposing reaction product enumeration."""

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .enumeration import enumerate_request
from .schemas import EnumerationRequest, EnumerationResponse
from .validation import RequestValidationError

app = FastAPI(
    title="Reaction Product Enumerator",
    description=(
        "Apply a known reaction SMARTS to reagent libraries and enumerate all "
        "valid product tuples with provenance. No yield prediction."
    ),
    version="1.0.0",
)


@app.exception_handler(RequestValidationError)
async def handle_validation_error(_: Request, exc: RequestValidationError):
    return JSONResponse(
        status_code=400,
        content={
            "error": "invalid input; whole batch rejected",
            "details": [item.model_dump() for item in exc.details],
        },
    )


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.post("/enumerate", response_model=EnumerationResponse)
async def enumerate_products(payload: EnumerationRequest) -> EnumerationResponse:
    return EnumerationResponse(**enumerate_request(payload))
