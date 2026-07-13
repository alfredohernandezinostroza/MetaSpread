"""Tests for the optional 3D spatial domain (space_dimensions == 3)."""
import numpy as np
import pytest

from metaspread import Config, run
from metaspread.cancermodel import CancerModel, _reflective_neighbor_sum
from metaspread.ndgrid import NDGrid


def small_3d(**overrides):
    base = Config.from_csv("simulations_configs.csv")
    cfg = base.copy(
        space_dimensions=3, gridsize=13, gridsize_z=13,
        number_of_initial_cells=10, n_center_points_for_tumor=10,
        n_center_points_for_Vessels=40, grids_number=2,
        extravasation_probs=[1.0], secondary_sites_vessels=[5],
        normal_vessels_primary=3, ruptured_vessels_primary=2,
    )
    return cfg.copy(**overrides) if overrides else cfg


def test_ndgrid_neighborhood_and_bounds():
    g = NDGrid((5, 5, 5))
    assert len(g.get_neighborhood((2, 2, 2), moore=False, include_center=False)) == 6
    with_center = g.get_neighborhood((2, 2, 2), moore=False, include_center=True)
    assert len(with_center) == 7 and (2, 2, 2) in with_center
    assert len(g.get_neighborhood((2, 2, 2), moore=True, include_center=False)) == 26
    # a corner has only 3 von-Neumann neighbours in bounds
    assert len(g.get_neighborhood((0, 0, 0), moore=False, include_center=False)) == 3
    assert g.out_of_bounds((5, 0, 0)) and not g.out_of_bounds((4, 4, 4))


def test_stencil_symmetry_3d():
    f = np.zeros((5, 5, 5))
    f[2, 2, 2] = 6.0
    s = _reflective_neighbor_sum(f)
    for nb in [(1, 2, 2), (3, 2, 2), (2, 1, 2), (2, 3, 2), (2, 2, 1), (2, 2, 3)]:
        assert s[nb] == 6.0
    assert s.sum() == 36.0


def test_3d_run_smoke():
    res = run(small_3d(), 5, 5, seed=1, save_path=None)
    m = res.model
    assert m.space_dimensions == 3
    assert m.mmp2[0].shape == (2, 13, 13, 13)
    cells = res.agent_data
    cells = cells[cells["Agent Type"] == "cell"]
    assert len(cells) > 0
    assert all(len(p) == 3 for p in cells["Position"])
    assert all(0 <= x < 13 and 0 <= y < 13 and 0 <= z < 13 for (x, y, z) in cells["Position"])


def test_3d_metastasis_occurs():
    # mobile cells + dense vessels + guaranteed survival so the cascade reaches a
    # secondary 3D site within a short run
    cfg = small_3d(
        dM=1e-3, dE=5e-4, vasculature_time=3,
        single_cell_survival=1.0, cluster_survival=1.0,
        number_of_initial_cells=30, n_center_points_for_tumor=30,
        n_center_points_for_Vessels=100, normal_vessels_primary=25,
        ruptured_vessels_primary=15, secondary_sites_vessels=[20],
    )
    res = run(cfg, 40, 40, seed=1, save_path=None)
    last = res.agent_data
    last = last[last["Step"] == last["Step"].max()]
    cells = last[last["Agent Type"] == "cell"]
    assert (cells["Grid"] == 2).sum() > 0


def test_3d_disk_save_and_reload(tmp_path):
    cfg = small_3d()
    run(cfg, 2, 2, seed=1, save_path=tmp_path)
    name = (f"Sim-max_steps-2-collection_period-2-"
            f"cells-{cfg.number_of_initial_cells}-grids_number-{cfg.grids_number}")
    sim_dir = tmp_path / "Simulations" / name
    # 3D fields are persisted as .npy
    assert len(list((sim_dir / "Mmp2").glob("*.npy"))) >= cfg.grids_number
    assert len(list((sim_dir / "Ecm").glob("*.npy"))) >= cfg.grids_number

    loaded = CancerModel(
        0, cfg.gridsize, cfg.gridsize, cfg.grids_number, 10, 10,
        sim_dir, loaded_simulation_path=sim_dir)
    assert loaded.space_dimensions == 3
    assert loaded.mmp2[0].shape == (2, 13, 13, 13)
    assert sum(1 for a in loaded.schedule.agents if a.agent_type == "cell") > 0
