"""Tests for the Phase 0 foundation: Config object, in-memory API, ensembles."""
import os

import pandas as pd
import pytest

import metaspread
from metaspread import Config, run, ensemble


def small_config():
    """A small, fast configuration for exercising the run machinery."""
    base = Config.from_csv("simulations_configs.csv")
    return base.copy(
        gridsize=41,
        number_of_initial_cells=10,
        n_center_points_for_tumor=10,
        grids_number=2,
        extravasation_probs=[1.0],
        secondary_sites_vessels=[10],
    )


def test_config_roundtrip(tmp_path):
    cfg = Config.from_csv("simulations_configs.csv")
    path = tmp_path / "roundtrip.csv"
    cfg.to_csv(path)
    reloaded = Config.from_csv(path)
    assert cfg.as_dict() == reloaded.as_dict()


def test_config_copy_is_independent():
    cfg = Config.from_csv("simulations_configs.csv")
    other = cfg.copy(carrying_capacity=99)
    assert other.carrying_capacity == 99
    assert cfg.carrying_capacity != 99  # original untouched


def test_in_memory_run_returns_data_and_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # any stray writes would land here
    metaspread.configs.generate_default_configs()
    cfg = small_config()
    results = run(cfg, max_steps=3, data_collection_period=3, seed=1, save_path=None)
    assert len(results.agent_data) > 0
    assert not os.path.exists(tmp_path / "Simulations")


def test_two_configs_same_process_are_independent():
    cfg = small_config()
    r1 = run(cfg.copy(number_of_initial_cells=10), 1, 1, seed=1, save_path=None)
    r2 = run(cfg.copy(number_of_initial_cells=25), 1, 1, seed=1, save_path=None)
    assert r1.model.number_of_initial_cells == 10
    assert r2.model.number_of_initial_cells == 25


def test_ensemble_reproducibility():
    cfg = small_config()
    same, _ = ensemble.run_ensemble(cfg, 3, 3, seeds=[5, 5])
    run0 = same[same["run_id"] == 0].reset_index(drop=True)
    run1 = same[same["run_id"] == 1].reset_index(drop=True)
    assert run0["Position"].tolist() == run1["Position"].tolist()

    diff, _ = ensemble.run_ensemble(cfg, 3, 3, seeds=[5, 6])
    d0 = diff[diff["run_id"] == 0].reset_index(drop=True)
    d1 = diff[diff["run_id"] == 1].reset_index(drop=True)
    assert d0["Position"].tolist() != d1["Position"].tolist()
