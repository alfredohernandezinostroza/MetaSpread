"""Tests for the Phase 2 organ-on-chip features (all off by default).

Feature 1: flow / advection of the diffusible fields (MMP2, oxygen).
"""
import numpy as np
import pytest

from metaspread import Config, run
from metaspread.cancermodel import _upwind_advection


def com_axis0(field):
    """Centre of mass of a 2D field along axis 0 (the x axis)."""
    weight = field.sum(axis=1)
    return float((weight * np.arange(field.shape[0])).sum() / field.sum())


# ---------------------------------------------------------------------------
# _upwind_advection: the operator itself
# ---------------------------------------------------------------------------
def test_upwind_advection_directionality():
    n = 6
    ramp = (np.arange(n).reshape(n, 1) * np.ones((1, n))).astype(float)  # field[x,y] = x

    # v > 0 uses a backward difference: derivative == slope (1) in the interior,
    # scaled by v; the upwind boundary row stays 0 (zero-gradient).
    a = _upwind_advection(ramp, [2.0, 0.0])
    assert np.allclose(a[1:, :], 2.0)
    assert np.allclose(a[0, :], 0.0)

    # v < 0 uses a forward difference: interior derivative 1*v, downwind edge 0.
    b = _upwind_advection(ramp, [-3.0, 0.0])
    assert np.allclose(b[:-1, :], -3.0)
    assert np.allclose(b[-1, :], 0.0)

    # a uniform field has zero gradient -> zero advection everywhere
    assert np.allclose(_upwind_advection(np.ones((n, n)), [5.0, -2.0]), 0.0)

    # zero velocity -> zero operator regardless of the field
    assert np.allclose(_upwind_advection(ramp, [0.0, 0.0]), 0.0)


def test_upwind_advection_is_3d_general():
    f = np.zeros((5, 5, 5))
    f[2, 2, 2] = 1.0
    out = _upwind_advection(f, [1.0, 0.0, 0.0])
    assert out.shape == f.shape
    # a positive x-velocity draws on the upwind (x-1) side of each cell
    assert out[2, 2, 2] == 1.0     # field[2]-field[1] = 1
    assert out[3, 2, 2] == -1.0    # field[3]-field[2] = -1


# ---------------------------------------------------------------------------
# integrated behaviour: a field blob is carried downstream
# ---------------------------------------------------------------------------
def test_advection_carries_blob_downstream():
    n = 41
    yy, xx = np.mgrid[0:n, 0:n]
    blob = np.exp(-(((xx - 20) ** 2 + (yy - 20) ** 2) / 8.0))
    tha, xha = 0.001, 0.005
    field = blob.copy()
    for _ in range(40):
        field = field - (tha / xha) * _upwind_advection(field, [0.6, 0.0])
    # first moment advects at ~v per unit time; downstream is +axis0
    assert com_axis0(field) > com_axis0(blob) + 1.0


# ---------------------------------------------------------------------------
# the wired path inside the model: flow shifts the MMP2 plume
# ---------------------------------------------------------------------------
def _flow_cfg(**overrides):
    base = Config.from_csv("simulations_configs.csv")
    cfg = base.copy(
        gridsize=61, grids_number=2, extravasation_probs=[1.0],
        secondary_sites_vessels=[5], number_of_initial_cells=20,
        n_center_points_for_tumor=20, n_center_points_for_Vessels=30,
    )
    return cfg.copy(**overrides) if overrides else cfg


def test_flow_advects_mmp2_field_in_solver():
    """The enable_flow branch inside calculate_environment moves the field.

    Driven deterministically: a Gaussian MMP2 blob is seeded and the real solver
    is stepped with the cell source zeroed, so only advection acts and the result
    can't be muddied by cell/RNG divergence between two stochastic runs.
    """
    cfg = _flow_cfg(enable_flow=True, flow_velocity=[0.6, 0.0, 0.0])
    model = run(cfg, 1, 1, seed=1, save_path=None).model

    # isolate advection: no recount, zero source everywhere
    model._recount_cells = lambda: None
    model.mesenchymal_count = [np.zeros_like(f[0]) for f in model.mmp2]
    model.epithelial_count = [np.zeros_like(f[0]) for f in model.mmp2]

    n = model.mmp2[0][0].shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    blob = np.exp(-(((xx - n // 2) ** 2 + (yy - n // 2) ** 2) / 6.0))
    for f in model.mmp2:
        f[0] = blob.copy(); f[1] = 0.0
    for e in model.ecm:
        e[0] = 0.5; e[1] = 0.5

    com0 = com_axis0(model.mmp2[0][0])
    for _ in range(25):
        model.calculate_environment(model.mmp2, model.ecm)
    com1 = com_axis0(model.mmp2[0][0])
    assert com1 > com0 + 0.5  # positive x-flow carried the plume to higher x


def test_flow_off_by_default():
    cfg = Config.from_csv("simulations_configs.csv")
    assert cfg.enable_flow is False
    # a default-config run must not touch the RNG / field path for flow: it just
    # has to complete without error and leave the flag inert.
    res = run(_flow_cfg(), 3, 3, seed=1, save_path=None)
    assert res.model.config.enable_flow is False
