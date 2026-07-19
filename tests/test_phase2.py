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


# ---------------------------------------------------------------------------
# Feature 2: shear-dependent survival of circulating clusters
# ---------------------------------------------------------------------------
from metaspread.cancermodel import get_cluster_survival_probability


def test_shear_scales_survival_probability():
    base = Config.from_csv("simulations_configs.csv")
    off = base.copy(enable_shear=False)
    on = base.copy(enable_shear=True, shear_stress=2.0, shear_death_coeff=0.5)
    single, cluster = (1, 0), (2, 1)

    # off: unchanged base probabilities
    assert get_cluster_survival_probability(single, off) == off.single_cell_survival
    assert get_cluster_survival_probability(cluster, off) == off.cluster_survival

    # on: scaled by exp(-coeff*stress) = exp(-1)
    factor = np.exp(-1.0)
    assert get_cluster_survival_probability(single, on) == pytest.approx(
        off.single_cell_survival * factor)
    assert get_cluster_survival_probability(cluster, on) == pytest.approx(
        off.cluster_survival * factor)

    # monotone: more shear -> less survival, and the factor stays in (0, 1]
    more = base.copy(enable_shear=True, shear_stress=5.0, shear_death_coeff=0.5)
    assert (get_cluster_survival_probability(single, more)
            < get_cluster_survival_probability(single, on))
    assert 0.0 < get_cluster_survival_probability(single, on) <= off.single_cell_survival


def test_shear_reduces_survivors_in_circulation():
    """Fewer circulating clusters survive under shear.

    Exercises the exact predicate the model uses to filter the vasculature
    (``random.random() < get_cluster_survival_probability(cluster, config)``) over
    a fixed RNG stream, so it is a fast, deterministic stand-in for the full
    metastatic cascade without the cost (and small-grid flooding) of running it.
    """
    import random

    off = Config.from_csv("simulations_configs.csv").copy(
        cluster_survival=0.8, enable_shear=False)
    on = off.copy(enable_shear=True, shear_stress=3.0, shear_death_coeff=0.5)
    clusters = [(2, 1)] * 500

    def survivors(cfg, seed):
        rng = random.Random(seed)
        return sum(1 for c in clusters
                   if rng.random() < get_cluster_survival_probability(c, cfg))

    # same seed, same clusters: shear must leave strictly fewer survivors
    assert survivors(on, 0) < survivors(off, 0)
    # and still let some through (factor exp(-1.5) ~ 0.22, base 0.8 -> ~0.18)
    assert survivors(on, 0) > 0


# ---------------------------------------------------------------------------
# Feature 3: device geometry (impassable walls)
# ---------------------------------------------------------------------------
from metaspread.configs import validate_configs
from metaspread.geometry import build_wall_mask


def test_wall_mask_parametric_channel_2d():
    base = Config.from_csv("simulations_configs.csv")
    cfg = base.copy(channel_axis=0, channel_margin=3, device_mask_path="")
    mask = build_wall_mask(cfg, (11, 11))
    # channel open along axis 0 (x): the first/last 3 columns (y) are walls,
    # the middle is fully open
    assert mask[:, :3].all() and mask[:, -3:].all()
    assert not mask[:, 3:-3].any()
    # margin 0 -> no walls at all (inert)
    assert not build_wall_mask(base.copy(channel_margin=0), (11, 11)).any()


def test_wall_mask_3d_and_mask_file(tmp_path):
    base = Config.from_csv("simulations_configs.csv")
    m = build_wall_mask(base.copy(channel_axis=2, channel_margin=1), (5, 5, 5))
    # channel open along z (axis 2): x and y edges are walls, a central z-column open
    assert m[:1].all() and m[-1:].all() and m[:, :1].all() and m[:, -1:].all()
    assert not m[2, 2, :].any()

    # a mask file overrides the parametric channel and is shape-checked
    arr = np.zeros((4, 4), dtype=int)
    arr[0, 0] = 1
    p = tmp_path / "mask.npy"
    np.save(p, arr)
    fm = build_wall_mask(base.copy(channel_margin=9, device_mask_path=str(p)), (4, 4))
    assert fm[0, 0] and fm.sum() == 1  # the file was used, not the channel
    with pytest.raises(ValueError):
        build_wall_mask(base.copy(device_mask_path=str(p)), (5, 5))  # shape mismatch


def test_geometry_off_by_default_and_validation():
    cfg = Config.from_csv("simulations_configs.csv")
    assert cfg.enable_device_geometry is False
    assert cfg.device_mask_path == ""  # empty cell round-trips to ""

    good = cfg.as_dict()
    good.update(enable_device_geometry=True, channel_axis=1, channel_margin=2)
    validate_configs(good)  # a valid channel config must not raise

    bad = cfg.as_dict()
    bad.update(enable_device_geometry=True, channel_axis=5)  # out of range for 2D
    with pytest.raises(ValueError):
        validate_configs(bad)


def _channel_cfg(**overrides):
    base = Config.from_csv("simulations_configs.csv")
    cfg = base.copy(
        gridsize=21, grids_number=2, extravasation_probs=[1.0],
        secondary_sites_vessels=[5], number_of_initial_cells=20,
        n_center_points_for_tumor=20, n_center_points_for_Vessels=40,
        enable_device_geometry=True, channel_axis=0, channel_margin=4,
    )
    return cfg.copy(**overrides) if overrides else cfg


def test_no_cancer_or_immune_agent_ever_enters_a_wall():
    # immune cells on too, so both movement paths are exercised against walls
    cfg = _channel_cfg(enable_immune=True, n_immune_cells=10, immune_kill_prob=0.0)
    res = run(cfg, 8, 1, seed=2, save_path=None)
    model = res.model

    # nothing alive sits on a wall (covers placement + movement end state)
    for a in model.schedule.agents:
        if a.agent_type in ("cell", "immune"):
            assert not model._is_wall(a.pos), f"{a.agent_type} on wall at {a.pos}"

    # and no cell/immune was ever recorded on a wall over the whole run
    hist = res.agent_data
    hist = hist[hist["Agent Type"].isin(["cell", "immune"])]
    assert len(hist) > 0
    assert not any(bool(model.wall_mask[tuple(pos)]) for pos in hist["Position"])
