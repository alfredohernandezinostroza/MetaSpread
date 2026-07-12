"""Tests for Phase 1: EMT/MET plasticity, oxygen field, immune agents.

All Phase 1 features are off by default; these tests enable them explicitly on a
small, fast configuration.
"""
import numpy as np
import pandas as pd
import pytest

import metaspread
from metaspread import Config, run
from metaspread.configs import PARAM_NAMES, DEFAULTS


def small_config(**overrides):
    base = Config.from_csv("simulations_configs.csv")
    cfg = base.copy(
        gridsize=41, number_of_initial_cells=10, n_center_points_for_tumor=10,
        grids_number=2, extravasation_probs=[1.0], secondary_sites_vessels=[10],
    )
    return cfg.copy(**overrides) if overrides else cfg


def _last_step_cells(results):
    df = results.agent_data
    last = df[df["Step"] == df["Step"].max()]
    return last[last["Agent Type"] == "cell"]


def test_emt_stochastic_converts_epithelial_to_mesenchymal():
    results = run(small_config(emt_prob=1.0, met_prob=0.0), 1, 1, seed=1, save_path=None)
    cells = _last_step_cells(results)
    assert (cells["Phenotype"] == "epithelial").sum() == 0
    assert (cells["Phenotype"] == "mesenchymal").sum() > 0


def test_met_stochastic_converts_mesenchymal_to_epithelial():
    results = run(small_config(emt_prob=0.0, met_prob=1.0), 1, 1, seed=1, save_path=None)
    cells = _last_step_cells(results)
    assert (cells["Phenotype"] == "mesenchymal").sum() == 0


def test_oxygen_field_forms_a_gradient():
    # start from an empty field so the vessel-driven gradient is unambiguous
    results = run(small_config(enable_oxygen=True, oxygen_initial=0.0),
                  20, 20, seed=1, save_path=None)
    oxygen = results.model.oxygen
    assert oxygen is not None
    primary = oxygen[0][0, :, :]
    assert primary.max() > primary.min()          # a gradient exists
    # oxygen is higher at a vessel than in the field on average
    vx, vy = results.model.grid_vessels_positions[0][0]
    assert primary[vx, vy] > primary.mean()


def test_immune_cells_reduce_cancer_burden():
    def final_cancer_count(cfg):
        return len(_last_step_cells(run(cfg, 12, 12, seed=1, save_path=None)))

    without_immune = final_cancer_count(small_config())
    # dense, fast-moving, lethal immune cells so contact with the tumour is certain
    with_immune = final_cancer_count(small_config(
        enable_immune=True, n_immune_cells=200,
        immune_kill_prob=1.0, immune_diff_coeff=0.02))
    assert with_immune < without_immune


def test_old_config_without_phase1_params_backfills(tmp_path):
    # a pre-Phase-1 configs.csv: core params + runtime rows, no Phase 1 params
    core = {name: DEFAULTS[name] for name in PARAM_NAMES if name not in (
        "emt_prob", "met_prob", "enable_hypoxia_emt", "hypoxia_threshold",
        "enable_oxygen", "d_oxygen", "oxygen_supply", "oxygen_consumption",
        "oxygen_initial", "oxygen_max", "enable_immune", "n_immune_cells",
        "immune_kill_prob", "immune_diff_coeff")}
    names = list(core.keys()) + ["max_steps", "data_collection_period"]
    values = list(core.values()) + [10, 10]
    path = tmp_path / "configs.csv"
    pd.DataFrame({"Names": names, "Values": values}).to_csv(path, index=False)

    with pytest.warns(UserWarning):
        cfg = Config.from_saved_simulation(path)
    assert cfg.enable_oxygen is False
    assert cfg.emt_prob == 0.0
    assert cfg.immune_diff_coeff == DEFAULTS["immune_diff_coeff"]
