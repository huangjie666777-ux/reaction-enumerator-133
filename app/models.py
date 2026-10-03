"""Pydantic request/response models."""
from __future__ import annotations

from pydantic import BaseModel, Field


MAX_COMBINATIONS = 500
MAX_PRODUCT_TUPLES_PER_COMBINATION = 100
MIN_REACTANT_POSITIONS = 1
MAX_REACTANT_POSITIONS = 3


class ReagentIn(BaseModel):
    id: str = Field(min_length=1, description="Unique reagent id within its reactant position")
    smiles: str = Field(min_length=1, description="Reagent structure as SMILES")


class EnumerationRequest(BaseModel):
    reaction_smarts: str = Field(min_length=1)
    reactants: list[list[ReagentIn]] = Field(description="Reagents grouped by reactant position")


class ProductOut(BaseModel):
    smiles: list[str]
    sources: list[list[str]]
    source_count: int


class FailureOut(BaseModel):
    reactant_ids: list[str]
    type: str
    message: str
    product_tuple_index: int | None = None


class EnumerationResponse(BaseModel):
    total_combinations: int
    processed_combinations: int
    no_match_combinations: int
    complete: bool
    products: list[ProductOut]
    failures: list[FailureOut]
    limits: dict[str, int]
