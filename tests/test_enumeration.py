from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

AMIDE = (
    "[C:1](=[O:2])([OH:3])[*:4].[NH2:5][*:6]"
    ">>[*:4][C:1](=[O:2])[N:5][*:6].[OH2:3]"
)


def post(payload):
    return client.post("/enumerate", json=payload)


def test_amide_cartesian_enumeration_and_sources():
    resp = post(
        {
            "reaction_smarts": AMIDE,
            "reagents": [
                [
                    {"id": "acetic", "smiles": "CC(=O)O"},
                    {"id": "benzoic", "smiles": "O=C(O)c1ccccc1"},
                ],
                [
                    {"id": "methylamine", "smiles": "CN"},
                    {"id": "ethylamine", "smiles": "CCN"},
                ],
            ],
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["total_combinations"] == 4
    assert body["processed_combinations"] == 4
    assert body["no_match_combinations"] == 0
    assert body["complete"] is True
    for product in body["products"]:
        assert len(product["product_tuple"]) == 2
        assert "O" in product["product_tuple"]
        assert len(product["sources"]) == 1
        assert len(product["sources"][0]) == 2
    keys = [tuple(p["product_tuple"]) for p in body["products"]]
    # Repeated requests yield identical ordering of products and sources.
    body2 = post(
        {
            "reaction_smarts": AMIDE,
            "reagents": [
                [
                    {"id": "acetic", "smiles": "CC(=O)O"},
                    {"id": "benzoic", "smiles": "O=C(O)c1ccccc1"},
                ],
                [
                    {"id": "methylamine", "smiles": "CN"},
                    {"id": "ethylamine", "smiles": "CCN"},
                ],
            ],
        }
    ).json()
    assert [tuple(p["product_tuple"]) for p in body2["products"]] == keys
    assert [p["sources"] for p in body2["products"]] == [
        p["sources"] for p in body["products"]
    ]


def test_equivalent_sites_dedup_within_combo():
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [[{"id": "cyclohexane", "smiles": "C1CCCCC1"}]],
        }
    )
    body = resp.json()
    assert body["complete"] is True
    assert len(body["products"]) == 1
    assert body["products"][0]["sources"] == [["cyclohexane"]]


def test_source_aggregation_across_combos():
    # Two reagent slots with the same identity template: identical SMILES with
    # different IDs funnel into one product; both sources must be listed once.
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [
                [
                    {"id": "c1", "smiles": "CCO"},
                    {"id": "c2", "smiles": "CCO"},
                ]
            ],
        }
    )
    body = resp.json()
    assert len(body["products"]) == 1
    assert body["products"][0]["sources"] == [["c1"], ["c2"]]


def test_no_match_counted_separately_from_invalid():
    resp = post(
        {
            "reaction_smarts": AMIDE,
            "reagents": [
                [{"id": "acetic", "smiles": "CC(=O)O"}],
                [{"id": "propane", "smiles": "CCC"}],
            ],
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["no_match_combinations"] == 1
    assert body["products"] == []
    assert body["failures"][0]["reason"] == "no_match"


def test_invalid_product_reported_per_combo():
    # Forcing a 5-coordinate carbon: product mol fails sanitization.
    resp = post(
        {
            "reaction_smarts": (
                "[Cl:1].[#6H4:2]"
                ">>[#6:2]([Cl:1])([Cl:1])([Cl:1])([Cl:1])([Cl:1])"
            ),
            "reagents": [
                [{"id": "cl", "smiles": "Cl"}],
                [{"id": "methane", "smiles": "C"}],
            ],
        }
    )
    body = resp.json()
    assert body["products"] == []
    reasons = {f["reason"] for f in body["failures"]}
    assert "invalid_product" in reasons
    detail = next(
        f["detail"] for f in body["failures"] if f["reason"] == "invalid_product"
    )
    assert detail


def test_stereoisomers_not_merged():
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [
                [
                    {"id": "r", "smiles": "C[C@H](O)C(=O)O"},
                    {"id": "s", "smiles": "C[C@@H](O)C(=O)O"},
                ]
            ],
        }
    )
    body = resp.json()
    assert len(body["products"]) == 2


def test_isotope_and_formal_charge_preserved():
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [[{"id": "ammonium15", "smiles": "[15NH4+]"}]],
        }
    )
    body = resp.json()
    smiles = body["products"][0]["product_tuple"][0]
    assert "15" in smiles
    assert "+" in smiles


def test_invalid_template_whole_batch_rejected_with_location():
    resp = post(
        {
            "reaction_smarts": "not_a_reaction",
            "reagents": [[{"id": "a", "smiles": "C"}]],
        }
    )
    assert resp.status_code == 400
    details = resp.json()["details"]
    assert all(d["location"] == "reaction_smarts" for d in details)


def test_invalid_reagent_smiles_points_to_position_and_index():
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [[{"id": "bad", "smiles": "notvalid(("}]],
        }
    )
    assert resp.status_code == 400
    d = resp.json()["details"][0]
    assert d["location"] == "reagent"
    assert d["position"] == 0
    assert d["index"] == 0
    assert d["reagent_id"] == "bad"


def test_position_count_mismatch_rejected():
    resp = post(
        {
            "reaction_smarts": AMIDE,
            "reagents": [[{"id": "a", "smiles": "CC(=O)O"}]],
        }
    )
    assert resp.status_code == 400
    assert "mismatch" in resp.json()["details"][0]["message"]


def test_duplicate_reagent_id_rejected():
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [
                [
                    {"id": "x", "smiles": "C"},
                    {"id": "x", "smiles": "O"},
                ]
            ],
        }
    )
    assert resp.status_code == 400
    assert "duplicate" in resp.json()["details"][0]["message"]


def test_over_500_combinations_pre_rejected():
    many = [{"id": f"r{i}", "smiles": "C"} for i in range(251)]
    resp = post(
        {
            "reaction_smarts": "[*:1].[*:2]>>[*:1][*:2]",
            "reagents": [many[:2], many],
        }
    )
    assert resp.status_code == 400
    assert "500" in resp.json()["details"][0]["message"]


def test_same_source_not_double_counted_and_input_unchanged():
    from rdkit import Chem

    smiles = "C1CCCCC1"
    before = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
    resp = post(
        {
            "reaction_smarts": "[*:1]>>[*:1]",
            "reagents": [[{"id": "cy", "smiles": smiles}]],
        }
    )
    body = resp.json()
    assert body["products"][0]["sources"].count(["cy"]) == 1
    after = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
    assert before == after


def test_tuple_limit_marked_incomplete():
    from rdkit import Chem

    from app.products import MAX_PRODUCT_TUPLES, run_combination
    from app.validation import ParsedReagent

    class FakeReaction:
        def RunReactants(self, mols, maxProducts=0):
            tuples = []
            for i in range(maxProducts):
                m = Chem.MolFromSmiles("[CH4]")
                m.GetAtomWithIdx(0).SetIsotope(i)
                tuples.append((m,))
            return tuples

    outcome = run_combination(
        FakeReaction(), (ParsedReagent("x", "C", Chem.MolFromSmiles("C")),)
    )
    assert outcome.truncated is True
    assert outcome.raw_tuple_count == MAX_PRODUCT_TUPLES + 1
    assert len(outcome.valid_keys) == MAX_PRODUCT_TUPLES


def test_health():
    assert client.get("/health").json()["status"] == "ok"
