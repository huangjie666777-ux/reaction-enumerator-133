from fastapi.testclient import TestClient

from app.main import app
from app.models import ReagentIn
from app.normalization import product_tuple_key

client = TestClient(app)

AMIDE_SMARTS = "[N:1].[C:2](=[O:3])Cl>>[N:1][C:2](=[O:3])"


def _reagent(rid, smiles):
    return {"id": rid, "smiles": smiles}


def test_basic_enumeration_and_sources():
    payload = {
        "reaction_smarts": AMIDE_SMARTS,
        "reactants": [
            [_reagent("a1", "CN"), _reagent("a2", "CCN")],
            [_reagent("b1", "CC(=O)Cl")],
        ],
    }
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_combinations"] == 2
    assert data["processed_combinations"] == 2
    assert data["no_match_combinations"] == 0
    assert data["complete"] is True
    smiles = {tuple(p["smiles"]) for p in data["products"]}
    assert smiles == {("CNC(C)=O",), ("CCNC(C)=O",)}
    by_smi = {p["smiles"][0]: p for p in data["products"]}
    assert by_smi["CNC(C)=O"]["sources"] == [["a1", "b1"]]
    assert by_smi["CNC(C)=O"]["source_count"] == 1
    assert by_smi["CCNC(C)=O"]["sources"] == [["a2", "b1"]]
    assert data["failures"] == []


def test_repeated_request_stable_order():
    payload = {
        "reaction_smarts": AMIDE_SMARTS,
        "reactants": [
            [_reagent("a1", "CN"), _reagent("a2", "CCN"), _reagent("a3", "N")],
            [_reagent("b1", "CC(=O)Cl"), _reagent("b2", "ClC(=O)C")],
        ],
    }
    first = client.post("/enumerate", json=payload).json()
    second = client.post("/enumerate", json=payload).json()
    assert first == second


def test_invalid_reaction_smarts_located():
    payload = {
        "reaction_smarts": "C>>(",
        "reactants": [[_reagent("a1", "C")]],
    }
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 400
    details = resp.json()["details"]
    assert any(d["location"] == "reaction_smarts" for d in details)


def test_position_count_mismatch_located():
    payload = {
        "reaction_smarts": AMIDE_SMARTS,
        "reactants": [[_reagent("a1", "CN")]],
    }
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 400
    assert "Expected 2" in resp.json()["details"][0]["message"]


def test_invalid_reagent_smiles_located_and_batch_rejected():
    payload = {
        "reaction_smarts": AMIDE_SMARTS,
        "reactants": [
            [_reagent("a1", "CN"), _reagent("a2", "not-a-smiles")],
            [_reagent("b1", "CC(=O)Cl")],
        ],
    }
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 400
    locs = {d["location"] for d in resp.json()["details"]}
    assert "reactants[0][1].smiles" in locs


def test_duplicate_reagent_id_rejected():
    payload = {
        "reaction_smarts": AMIDE_SMARTS,
        "reactants": [
            [_reagent("dup", "CN"), _reagent("dup", "CCN")],
            [_reagent("b1", "CC(=O)Cl")],
        ],
    }
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 400
    assert any("Duplicate" in d["message"] for d in resp.json()["details"])


def test_combination_limit_pre_rejected():
    amines = [_reagent(f"a{i}", "CN") for i in range(26)]
    acids = [_reagent(f"b{i}", "CC(=O)Cl") for i in range(20)]
    payload = {"reaction_smarts": AMIDE_SMARTS, "reactants": [amines, acids]}
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 400
    msg = resp.json()["details"][0]["message"]
    assert "520" in msg and "500" in msg


def test_no_match_counted():
    payload = {
        "reaction_smarts": AMIDE_SMARTS,
        "reactants": [
            [_reagent("a1", "c1ccccc1")],
            [_reagent("b1", "CC(=O)Cl")],
        ],
    }
    data = client.post("/enumerate", json=payload).json()
    assert data["total_combinations"] == 1
    assert data["no_match_combinations"] == 1
    assert data["products"] == []


def test_stereoisomers_not_merged():
    payload = {
        "reaction_smarts": "[C:1]>>[C:1]",
        "reactants": [
            [
                _reagent("s1", "N[C@@H](C)C(=O)O"),
                _reagent("s2", "N[C@H](C)C(=O)O"),
            ]
        ],
    }
    data = client.post("/enumerate", json=payload).json()
    smiles = {p["smiles"][0] for p in data["products"]}
    assert any("@" in s for s in smiles)
    assert len(smiles) == 2


def test_isotope_and_charge_preserved_and_maps_removed():
    payload = {
        "reaction_smarts": "[N:1]>>[N:1]",
        "reactants": [[_reagent("i1", "[15NH4+]")]],
    }
    data = client.post("/enumerate", json=payload).json()
    smi = data["products"][0]["smiles"][0]
    assert "15" in smi and "+" in smi


def test_multiple_product_templates_keep_full_tuple():
    payload = {
        "reaction_smarts": "[C:1].[O:2]>>[C:1].[O:2]",
        "reactants": [[_reagent("c1", "CC")], [_reagent("o1", "CO")]],
    }
    data = client.post("/enumerate", json=payload).json()
    assert data["failures"] == []
    for product in data["products"]:
        assert len(product["smiles"]) == 2


def test_three_reactants():
    payload = {
        "reaction_smarts": "[C:1].[C:2].[C:3]>>[C:1][C:2][C:3]",
        "reactants": [
            [_reagent("x", "C")],
            [_reagent("y", "C")],
            [_reagent("z", "C")],
        ],
    }
    data = client.post("/enumerate", json=payload).json()
    assert data["total_combinations"] == 1
    assert len(data["products"][0]["smiles"]) == 1


def test_invalid_product_recorded_not_mixed():
    payload = {
        "reaction_smarts": "[C:1]>>[C:1][Cl][Cl]",
        "reactants": [[_reagent("x", "C")]],
    }
    data = client.post("/enumerate", json=payload).json()
    assert data["products"] == []
    assert data["failures"][0]["type"] == "invalid_product"
    assert data["failures"][0]["reactant_ids"] == ["x"]


def test_shared_product_aggregates_unique_sources():
    payload = {
        "reaction_smarts": "[C:1]>>[C:1]",
        "reactants": [[_reagent("a", "CC"), _reagent("b", "CC")]],
    }
    data = client.post("/enumerate", json=payload).json()
    sources = data["products"][0]["sources"]
    assert sources == [["a"], ["b"]]
    assert data["products"][0]["source_count"] == 2


def test_too_many_tuples_marks_incomplete(monkeypatch):
    from app import reaction_runner

    monkeypatch.setattr(reaction_runner, "MAX_PRODUCT_TUPLES_PER_COMBINATION", 1)
    payload = {
        "reaction_smarts": "[C:1]>>[C:1]",
        "reactants": [[_reagent("p", "CCC")]],
    }
    data = client.post("/enumerate", json=payload).json()
    assert data["complete"] is False
    assert any(f["type"] == "too_many_product_tuples" for f in data["failures"])


def test_product_tuple_key_order_independent_with_multiplicity():
    assert product_tuple_key(["B", "A"]) == product_tuple_key(["A", "B"])
    assert product_tuple_key(["A", "A", "B"]) != product_tuple_key(["A", "B"])
    assert product_tuple_key(["A@", "A@@"]) != product_tuple_key(["A@", "A@"])


def test_input_molecules_not_mutated():
    from rdkit import Chem

    from app.aggregation import aggregate
    from app.validation import validate_request

    original = Chem.MolFromSmiles("CN")
    before = Chem.MolToSmiles(original)
    rxn, positions = validate_request(
        AMIDE_SMARTS,
        [
            [ReagentIn(id="a1", smiles="CN")],
            [ReagentIn(id="b1", smiles="CC(=O)Cl")],
        ],
    )
    aggregate(rxn, positions)
    assert Chem.MolToSmiles(original) == before


def test_empty_position_rejected():
    payload = {"reaction_smarts": "[C:1]>>[C:1]", "reactants": [[]]}
    resp = client.post("/enumerate", json=payload)
    assert resp.status_code == 400
