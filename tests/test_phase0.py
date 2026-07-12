"""Tests for the Phase 0 foundation: Config object, in-memory API, ensembles."""
import os

import pandas as pd
import pytest

import metaspread
import metaspread.configs
from metaspread import Config, run, ensemble
from metaspread.configs import PARAM_NAMES, validate_configs


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


def test_disk_backed_run_writes_artifacts(tmp_path):
    cfg = small_config()
    results = run(cfg, max_steps=2, data_collection_period=2, seed=1, save_path=tmp_path)
    name = (f"Sim-max_steps-2-collection_period-2-"
            f"cells-{cfg.number_of_initial_cells}-grids_number-{cfg.grids_number}")
    sim_dir = tmp_path / "Simulations" / name
    assert (sim_dir / "CellsData.csv").is_file()
    for sub in ["Mmp2", "Ecm", "Vasculature", "Time when grids were populated"]:
        assert (sim_dir / sub).is_dir()
    # saved configs.csv carries the core params plus the two runtime rows
    saved = pd.read_csv(sim_dir / "configs.csv")
    assert len(saved) == len(PARAM_NAMES) + 2
    assert len(results.agent_data) > 0


def test_from_saved_simulation_roundtrip(tmp_path):
    cfg = small_config()
    run(cfg, max_steps=2, data_collection_period=2, seed=1, save_path=tmp_path)
    name = (f"Sim-max_steps-2-collection_period-2-"
            f"cells-{cfg.number_of_initial_cells}-grids_number-{cfg.grids_number}")
    saved_cfg = Config.from_saved_simulation(tmp_path / "Simulations" / name / "configs.csv")
    # core params preserved; runtime rows dropped
    assert saved_cfg.as_dict() == cfg.as_dict()
    assert not hasattr(saved_cfg, "max_steps")


@pytest.mark.parametrize("overrides", [
    {"extravasation_probs": [0.5]},                       # does not sum to 1
    {"extravasation_probs": [0.5, 0.5]},                  # wrong length vs grids_number
    {"secondary_sites_vessels": [10, 10]},                # wrong length vs grids_number
    {"mesenchymal_proportion": 0.5, "epithelial_proportion": 0.4},  # proportions != 1
    {"number_of_initial_cells": 41},                      # > n_center_points * capacity (10*4)
])
def test_validate_configs_rejects_bad_values(overrides):
    values = small_config().as_dict()
    values.update(overrides)
    with pytest.raises(ValueError):
        validate_configs(values)


def test_ensemble_param_overrides():
    cfg = small_config()
    baseline_capacity = cfg.carrying_capacity
    combined, results = ensemble.run_ensemble(
        cfg, 1, 1, seeds=[1],
        param_overrides=[{"carrying_capacity": 4}, {"carrying_capacity": 8}],
    )
    assert sorted(combined["run_id"].unique().tolist()) == [0, 1]
    assert set(combined["carrying_capacity"].unique()) == {"4", "8"}
    assert results[0].config.carrying_capacity == 4
    assert results[1].config.carrying_capacity == 8
    # base config untouched by the sweep
    assert cfg.carrying_capacity == baseline_capacity


def test_determinism_same_seed_identical():
    cfg = small_config()
    r1 = run(cfg, 2, 2, seed=3, save_path=None)
    r2 = run(cfg, 2, 2, seed=3, save_path=None)
    assert r1.agent_data["Position"].tolist() == r2.agent_data["Position"].tolist()


def test_publish_to_module_contract():
    cfg = Config.from_csv("simulations_configs.csv")
    cfg.publish_to_module()
    # the backward-compat contract datagenerator/graphgenerator rely on
    assert metaspread.configs.gridsize == cfg.gridsize
    assert metaspread.configs.carrying_capacity == cfg.carrying_capacity
