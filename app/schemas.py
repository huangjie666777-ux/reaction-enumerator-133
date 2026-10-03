"""Pydantic request/response models for the enumeration API."""

from typing import List

from pydantic import BaseModel, Field


class Reagent(BaseModel):
    id: str = Field(..., min_length=1, description="Unique reagent ID within its slot")
    smiles: str = Field(..., min_length=1, description="Reagent SMILES")


class EnumerationRequest(BaseModel):
    reaction_smarts: str = Field(..., min_length=1)
    # Outer list = reactant positions (1..3), inner list = reagents at that position.
    reagents: List[List[Reagent]]


class ProductRecord(BaseModel):
    # Canonical isomeric SMILES of every component of one product tuple,
    # sorted (order-independent, multiplicity preserved).
    product_tuple: List[str]
    # Every source reagent-ID combination producing this tuple, first-seen order.
    sources: List[List[str]]


class FailureRecord(BaseModel):
    reagent_ids: List[str]
    reason: str  # "no_match" or "invalid_product"
    detail: str = ""
    truncated_product_tuples: bool = False


class EnumerationResponse(BaseModel):
    total_combinations: int
    processed_combinations: int
    no_match_combinations: int
    complete: bool
    products: List[ProductRecord]
    failures: List[FailureRecord]


class ValidationErrorItem(BaseModel):
    location: str  # "reaction_smarts" or "reagent" or "request"
    position: int = -1
    index: int = -1
    reagent_id: str = ""
    message: str


class ValidationErrorResponse(BaseModel):
    error: str
    details: List[ValidationErrorItem]
