from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

import plasrip
from plasrip import api, core, pipeline
from plasrip.schemas import CORE_SUMMARY_COLUMNS, SUMMARY_COLUMNS


def _fake_db(tmp_path: Path):
    return SimpleNamespace(
        root=tmp_path / "db",
        database_version="v2.0",
        schema_version="2.0",
        ptu_metadata=pd.DataFrame([{"PTU_ID": "PTU000001"}]),
        pgg_metadata=pd.DataFrame([{"PGG_ID": "PGG-000001"}]),
    )


def test_step7_public_api_exports_database_helpers():
    assert callable(plasrip.load_reference_database)
    assert callable(plasrip.validate_reference_database)
    assert plasrip.PlasRipReferenceDatabase is api.PlasRipReferenceDatabase


def test_step7_api_database_wrappers_delegate(monkeypatch):
    seen = {}

    def fake_load(database_dir, check_sketches=True):
        seen["load"] = (database_dir, check_sketches)
        return "DB"

    def fake_validate(database_dir, check_sketches=False, check_all_sketches=False):
        seen["validate"] = (database_dir, check_sketches, check_all_sketches)
        return {"ok": True}

    monkeypatch.setattr(api, "load_database", fake_load)
    monkeypatch.setattr(api, "validate_database_dir", fake_validate)

    assert api.load_reference_database("/db", check_sketches=False) == "DB"
    assert api.validate_reference_database("/db", check_sketches=True, check_all_sketches=True) == {"ok": True}
    assert seen["load"] == ("/db", False)
    assert seen["validate"] == ("/db", True, True)


def test_step7_core_accepts_preloaded_database_without_reloading(monkeypatch, tmp_path):
    db = _fake_db(tmp_path)
    assignment_indexes = object()
    ptu_lookup = {"PTU000001": {}}
    ptu_annotation = pd.DataFrame([{"PTU": "PTU000001"}])

    def should_not_load(*args, **kwargs):
        raise AssertionError("load_database must not be called when _database is supplied")

    monkeypatch.setattr(core, "load_database", should_not_load)

    def fake_predict(**kwargs):
        return {
            "hostrange": pd.DataFrame([{
                "PLASMID": "query1",
                "Predicted_host_top1": "Escherichia coli",
            }]),
            "top_hits": pd.DataFrame(),
            "taxonomy": pd.DataFrame(),
            "set_df": pd.DataFrame(),
        }

    seen = {}

    def fake_assign(**kwargs):
        seen["database"] = kwargs.get("_database")
        seen["indexes"] = kwargs.get("_assignment_indexes")
        seen["lookup"] = kwargs.get("_ptu_meta_lookup")
        return pd.DataFrame([{
            "query": "query1",
            "ptu_assigned": "PTU000001",
            "PGG_ID": "PGG-000001",
            "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
            "n_ptu_hits": 1,
            "how_assigned": "closest_rep",
            "closest_rep_dist": 0.01,
        }])

    monkeypatch.setattr(core, "run_host_prediction", fake_predict)
    monkeypatch.setattr(core, "run_ptu_assignment", fake_assign)

    result = core.run_core_analysis(
        input_fastas=["query1.fasta"],
        model_dir=None,
        database_dir="/ignored",
        out_dir=str(tmp_path / "core"),
        quiet=True,
        overwrite=True,
        _database=db,
        _assignment_indexes=assignment_indexes,
        _ptu_meta_lookup=ptu_lookup,
        _ptu_annotation=ptu_annotation,
    )

    assert seen["database"] is db
    assert seen["indexes"] is assignment_indexes
    assert seen["lookup"] is ptu_lookup
    assert list(result["core_summary"].columns) == CORE_SUMMARY_COLUMNS
    assert result["core_summary"].iloc[0]["PTU_assignment_method"] == "closest_rep"


def test_step7_pipeline_loads_database_once_and_reuses_it(monkeypatch, tmp_path):
    db = _fake_db(tmp_path)
    calls = {"load": 0, "core": []}
    indexes = object()
    lookup = object()
    annotation = pd.DataFrame([{"PTU": "PTU000001"}])

    def fake_load(database_dir, check_sketches=True):
        calls["load"] += 1
        assert database_dir == "/db"
        assert check_sketches is True
        return db

    monkeypatch.setattr(pipeline, "load_database", fake_load)
    monkeypatch.setattr(pipeline, "build_assignment_indexes", lambda d: indexes)
    monkeypatch.setattr(pipeline, "ptu_metadata_lookup", lambda d: lookup)
    monkeypatch.setattr(pipeline, "_ptu_annotation_table", lambda d: annotation)

    def fake_core(**kwargs):
        calls["core"].append(kwargs)
        plasmid = Path(kwargs["input_fastas"][0]).stem
        return {
            "hostrange": pd.DataFrame(),
            "top_hits": pd.DataFrame(),
            "taxonomy": pd.DataFrame(),
            "set_df": pd.DataFrame(),
            "ptu": pd.DataFrame(),
            "ptu_hits": pd.DataFrame(),
            "core_summary": pd.DataFrame([{
                "PLASMID": plasmid,
                "PTU": "PTU000001",
                "PGG_ID": "PGG-000001",
                "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
                "PTU_assignment_method": "closest_rep",
                "PTU_candidate_hits": 1,
                "PTU_closest_rep_dist": 0.01,
            }]),
        }

    def fake_annotation(**kwargs):
        plasmid = Path(kwargs["input_fastas"][0]).stem
        return {
            "amr": pd.DataFrame(),
            "replicon": pd.DataFrame(),
            "mobility": pd.DataFrame(),
            "annotation_summary": pd.DataFrame([{
                "PLASMID": plasmid,
                "AMR_hits": 0,
                "AMR_genes": "",
                "REPLICON_hits": 0,
                "REPLICON_genes": "",
                "MOB": "NT",
            }]),
        }

    monkeypatch.setattr(pipeline, "run_core_analysis", fake_core)
    monkeypatch.setattr(pipeline, "run_annotation", fake_annotation)

    result = pipeline.run_pipeline(
        input_fastas=["a.fasta", "b.fasta"],
        out_dir=str(tmp_path / "run"),
        model_dir=None,
        database_dir="/db",
        quiet=True,
        overwrite=True,
    )

    assert calls["load"] == 1
    assert len(calls["core"]) == 2
    for kwargs in calls["core"]:
        assert kwargs["_database"] is db
        assert kwargs["_assignment_indexes"] is indexes
        assert kwargs["_ptu_meta_lookup"] is lookup
        assert kwargs["_ptu_annotation"] is annotation
        assert kwargs["database_dir"] == str(db.root)

    assert len(result["summary"]) == 2
    assert list(result["summary"].columns) == SUMMARY_COLUMNS
    assert set(result["summary"]["PLASMID"]) == {"a", "b"}
    assert set(result["summary"]["PGG_ID"]) == {"PGG-000001"}
