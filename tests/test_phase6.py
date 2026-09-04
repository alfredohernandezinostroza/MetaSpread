"""Tests for Phase 6: the presentation-grade ("showcase") viewer.

Like Phase 5 this is postprocessing only — read-only over saved artifacts — so
it cannot affect the byte-identity gate. What these pin down is that the pretty
view still shows *real* data: one animation frame per saved step, a device frame
rebuilt from the run's own config, and an export that is genuinely standalone.

plotly lives in the optional 'viz' extra, so the whole module skips on a base
install; run it with `pixi run -e viz pytest`.
"""
import re

import matplotlib
matplotlib.use("Agg")  # headless: never try to open a window during tests

import pytest

pytest.importorskip("plotly", reason="the showcase view needs the optional 'viz' extra")

from metaspread import Config, run, viewer3d


def _small_config(**overrides):
    base = Config.from_csv("simulations_configs.csv")
    cfg = base.copy(
        gridsize=41, number_of_initial_cells=10, n_center_points_for_tumor=10,
        grids_number=2, extravasation_probs=[1.0], secondary_sites_vessels=[10],
    )
    return cfg.copy(**overrides) if overrides else cfg


def _sim_name(cfg, max_steps, period):
    return (f"Sim-max_steps-{max_steps}-collection_period-{period}-"
            f"cells-{cfg.number_of_initial_cells}-grids_number-{cfg.grids_number}")


def _run_3d(tmp_path, max_steps=4, period=2):
    cfg = _small_config(space_dimensions=3, gridsize=11, gridsize_z=5,
                        number_of_initial_cells=8, n_center_points_for_tumor=8,
                        enable_oxygen=True)
    run(cfg, max_steps, period, seed=1, save_path=tmp_path)
    return str(tmp_path / "Simulations" / _sim_name(cfg, max_steps, period))


def test_showcase_figure_animates_every_saved_step(tmp_path):
    sim = _run_3d(tmp_path)
    data = viewer3d.load_simulation(sim)
    fig = viewer3d.showcase_figure(data, grid=1)

    # one frame per saved step, in order, addressed by step number
    assert [f.name for f in fig.frames] == [str(s) for s in data["steps"]]
    # frames replace traces by index, so every frame must carry the same count
    counts = {len(f.data) for f in fig.frames}
    assert len(counts) == 1, "frames disagree on trace count; animation would corrupt"
    # ...and the base figure starts with exactly those animated traces
    assert len(fig.data) >= counts.pop()
    # the slider exposes each step
    assert [s.label for s in fig.layout.sliders[0].steps] == [str(s) for s in data["steps"]]


def test_showcase_plots_only_real_agents(tmp_path):
    # the glow halo is a second copy of the same points - it must not invent any
    sim = _run_3d(tmp_path)
    data = viewer3d.load_simulation(sim)
    step = data["steps"][-1]
    groups = viewer3d.agent_groups(data["cells"], step, 1)

    glowing = viewer3d._agent_traces(groups, glow=True)
    plain = viewer3d._agent_traces(groups, glow=False)
    assert len(glowing) == 2 * len(plain)          # halo + core per group
    for trace in glowing:
        label = trace.name
        assert list(trace.x) == list(groups[label][0])


def test_device_frame_is_rebuilt_from_the_saved_config():
    # device geometry is not persisted; the frame comes from build_wall_mask,
    # a pure function of the config. No simulation run needed.
    cfg = _small_config()

    off = viewer3d.device_frame_traces(cfg.copy(enable_device_geometry=False), (21, 21))
    assert [t.name for t in off] == ["Domain"]

    on = viewer3d.device_frame_traces(
        cfg.copy(enable_device_geometry=True, channel_axis=0, channel_margin=4), (21, 21))
    assert [t.name for t in on] == ["Domain", "Device channel"]
    # the channel box spans the open corridor: 4 cells of wall on each side of y
    ys = [y for y in on[1].y if y is not None]
    assert (min(ys), max(ys)) == (4, 16)


def test_device_frame_ignores_geometry_with_no_walls():
    # enable_device_geometry with channel_margin 0 is an open domain: nothing to draw
    cfg = _small_config(enable_device_geometry=True, channel_margin=0)
    assert [t.name for t in viewer3d.device_frame_traces(cfg, (21, 21))] == ["Domain"]


def test_agent_counts_match_the_plotted_groups(tmp_path):
    sim = _run_3d(tmp_path)
    data = viewer3d.load_simulation(sim)
    step = data["steps"][-1]

    counts = viewer3d.agent_counts(data["cells"], step, 1)
    groups = viewer3d.agent_groups(data["cells"], step, 1)
    for label, (xs, _, _) in groups.items():
        assert counts[label] == len(xs)
    assert counts["Cancer cells"] == counts["Mesenchymal"] + counts["Epithelial"]
    assert counts["Cancer cells"] > 0


def test_showcase_works_for_2d_runs(tmp_path):
    # 2D positions have no z; the showcase view must still build (z collapses to 0)
    cfg = _small_config()
    run(cfg, 4, 2, seed=1, save_path=tmp_path)
    sim = str(tmp_path / "Simulations" / _sim_name(cfg, 4, 2))

    data = viewer3d.load_simulation(sim)
    assert data["space_dimensions"] == 2
    fig = viewer3d.showcase_figure(data, grid=1)
    assert len(fig.frames) == len(data["steps"])
    zs = [z for trace in fig.data for z in (trace.z or []) if z is not None]
    assert zs and set(zs) == {0}


def test_export_showcase_html_is_self_contained(tmp_path):
    # the point of the export is that it opens offline as one file
    sim = _run_3d(tmp_path)
    out = tmp_path / "showcase.html"
    viewer3d.export_showcase_html(sim, str(out))

    html = out.read_text(encoding="utf-8")
    assert out.stat().st_size > 500_000, "plotly.js does not look inlined"
    # what actually matters: the page loads nothing over the network. (The
    # bundle still *mentions* cdn.plot.ly as its topojson default, which only
    # geo/choropleth traces would ever fetch - we render Scatter3d.)
    assert not re.findall(r"<script[^>]*\ssrc=", html), "remote script tag"
    assert not re.findall(r"<link[^>]*\shref=", html), "remote stylesheet"
    assert "scatter3d" in html.lower()
