"""Data + figure helpers for the interactive 3D viewer.

The marimo UI (``viewer3d_app.py``) is a thin wrapper over these functions; all
the substantive logic lives here so it can be imported and unit-tested without
the optional ``viz`` dependencies (marimo, plotly). plotly is imported lazily
inside the figure builders, so the data helpers work with the base install.

Local use::

    pixi run -e viz marimo run metaspread/viewer3d_app.py -- --sim Simulations/<name>
"""
import ast
import os

import numpy as np
import pandas as pd

import metaspread.configs

# label -> (plotly colour, marker symbol) for the agent point cloud
AGENT_STYLE = {
    "Mesenchymal": ("blue", "circle"),
    "Epithelial": ("orange", "diamond"),
    "Vessel": ("red", "circle"),
    "Ruptured vessel": ("darkred", "x"),
    "Immune": ("green", "cross"),
}

# Same agent semantics, retuned for the dark "showcase" scene: the plain
# palette (and darkred in particular) washes out against a dark background.
SHOWCASE_STYLE = {
    "Mesenchymal": ("#46b8f4", "circle"),
    "Epithelial": ("#f7a63e", "diamond"),
    "Vessel": ("#ff5c7a", "circle"),
    "Ruptured vessel": ("#ff2d2d", "x"),
    "Immune": ("#2dd4c4", "cross"),
}

THEMES = {"plain": AGENT_STYLE, "showcase": SHOWCASE_STYLE}

# Fixed trace order. Animation frames address traces by index, so every frame
# must emit the same groups in the same order — even the empty ones.
AGENT_LABELS = tuple(AGENT_STYLE)

_DARK_BG = "#090d14"


def load_simulation(simulation_path):
    """Load a saved simulation and return the metadata the viewer needs.

    Publishes the run's configs onto ``metaspread.configs`` (like the other
    post-processors) and reads ``CellsData.csv`` once.
    """
    configs_path = os.path.join(simulation_path, "configs.csv")
    metaspread.configs.load_simulation_configs_for_data_generation(configs_path)
    cells = pd.read_csv(os.path.join(simulation_path, "CellsData.csv"),
                        converters={"Position": ast.literal_eval})
    return {
        "path": simulation_path,
        "cells": cells,
        "steps": sorted(cells["Step"].unique().tolist()),
        "grids": sorted(cells["Grid"].unique().tolist()),
        "gridsize": metaspread.configs.gridsize,
        "gridsize_z": int(getattr(metaspread.configs, "gridsize_z", 1)),
        "space_dimensions": int(getattr(metaspread.configs, "space_dimensions", 2)),
        "fields": [f for f in ("Oxygen", "Mmp2", "Ecm")
                   if os.path.isdir(os.path.join(simulation_path, f))],
    }


def agent_groups(cells, step, grid):
    """Split the agents at ``(step, grid)`` into plottable groups.

    Returns ``{label: (xs, ys, zs)}``; ``zs`` is 0 for 2D positions so the same
    3D scatter renders both 2D and 3D runs.
    """
    df = cells[(cells["Step"] == step) & (cells["Grid"] == grid)]

    def xyz(sub):
        pos = list(sub["Position"])
        xs = [p[0] for p in pos]
        ys = [p[1] for p in pos]
        zs = [p[2] if len(p) > 2 else 0 for p in pos]
        return xs, ys, zs

    return {
        "Mesenchymal": xyz(df[df["Phenotype"] == "mesenchymal"]),
        "Epithelial": xyz(df[df["Phenotype"] == "epithelial"]),
        "Vessel": xyz(df[(df["Agent Type"] == "vessel") & (df["Ruptured"] == False)]),
        "Ruptured vessel": xyz(df[(df["Agent Type"] == "vessel") & (df["Ruptured"] == True)]),
        "Immune": xyz(df[df["Agent Type"] == "immune"]),
    }


def load_field(simulation_path, field, grid, step):
    """Load a saved field for ``(grid, step)`` as an ndarray (2D CSV or 3D .npy)."""
    directory = os.path.join(simulation_path, field)
    npy = os.path.join(directory, f"{field}-{grid}grid-{step}step.npy")
    if os.path.isfile(npy):
        return np.load(npy)
    csv = os.path.join(directory, f"{field}-{grid}grid-{step}step.csv")
    return pd.read_csv(csv, index_col=0).values


def field_slice(simulation_path, field, grid, step, z):
    """Return one 2D z-slice of a field (or the whole 2D field if it is not 3D)."""
    arr = load_field(simulation_path, field, grid, step)
    if arr.ndim == 3:
        z = int(np.clip(z, 0, arr.shape[2] - 1))
        return arr[:, :, z]
    return arr


def agent_scatter_figure(groups):
    """Build a plotly 3D scatter (orbit/zoom) of the agent groups."""
    import plotly.graph_objects as go

    fig = go.Figure()
    for label, (xs, ys, zs) in groups.items():
        if not xs:
            continue
        colour, symbol = AGENT_STYLE[label]
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="markers", name=label,
            marker=dict(size=3, color=colour, symbol=symbol, opacity=0.8)))
    fig.update_layout(scene=dict(xaxis_title="x", yaxis_title="y", zaxis_title="z",
                                 aspectmode="data"),
                      height=650, legend=dict(itemsizing="constant"))
    return fig


def field_heatmap_figure(slice2d, title=""):
    """Build a plotly heatmap of a single field slice."""
    import plotly.graph_objects as go

    fig = go.Figure(go.Heatmap(z=np.asarray(slice2d).T, colorscale="Viridis"))
    fig.update_layout(title=title, height=420,
                      yaxis=dict(scaleanchor="x", scaleratio=1))
    return fig


# --------------------------------------------------------------------------
# Showcase view: the same real data, rendered for presentation (Phase 6).
#
# Everything below plots only what the simulation actually saved. The two
# concessions to looks are the glow halo (a second, larger, translucent copy of
# each marker) and the dark scene; neither invents an agent or a position.
# --------------------------------------------------------------------------


def agent_counts(cells, step, grid):
    """Per-group agent counts at ``(step, grid)`` — the numbers behind the HUD."""
    counts = {label: len(xs) for label, (xs, _, _) in agent_groups(cells, step, grid).items()}
    counts["Cancer cells"] = counts["Mesenchymal"] + counts["Epithelial"]
    return counts


def _box_edges(lo, hi):
    """xs/ys/zs tracing the 12 edges of an axis-aligned box as one gapped polyline."""
    corners = [(x, y, z) for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]
    xs, ys, zs = [], [], []
    for a in range(len(corners)):
        for b in range(a + 1, len(corners)):
            # an edge joins corners differing in exactly one coordinate
            if sum(ca != cb for ca, cb in zip(corners[a], corners[b])) != 1:
                continue
            for axis, seq in enumerate((xs, ys, zs)):
                seq.extend([corners[a][axis], corners[b][axis], None])
    return xs, ys, zs


def _spatial_shape(data):
    """The field shape (W, H) or (W, H, D) implied by a loaded simulation."""
    shape = (data["gridsize"], data["gridsize"])
    if data["space_dimensions"] == 3:
        shape = shape + (data["gridsize_z"],)
    return shape


def device_frame_traces(config, shape, theme="showcase"):
    """Wireframe traces for the simulation domain, plus the device channel.

    The wall mask is *rebuilt* from the saved config — ``build_wall_mask`` is a
    pure function of it — so nothing extra has to have been persisted by the run.
    When device geometry is off, only the domain box is drawn.

    The channel is drawn as the bounding box of the open (non-wall) region, which
    works for both the parametric channel and an arbitrary ``device_mask_path``
    mask. Rendering the full wall volume is a documented follow-up.
    """
    import plotly.graph_objects as go

    shape = tuple(int(s) for s in shape)
    dark = theme == "showcase"
    hi = (shape[0] - 1, shape[1] - 1, (shape[2] - 1) if len(shape) > 2 else 0)

    xs, ys, zs = _box_edges((0, 0, 0), hi)
    traces = [go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", name="Domain", hoverinfo="skip",
        line=dict(color="#7896b4" if dark else "#b0b8c0", width=1))]

    if not getattr(config, "enable_device_geometry", False):
        return traces

    from metaspread.geometry import build_wall_mask
    mask = build_wall_mask(config, shape)
    open_idx = np.nonzero(~mask)
    if not mask.any() or len(open_idx[0]) == 0:
        return traces          # geometry on but no walls: nothing extra to show

    olo = [int(a.min()) for a in open_idx]
    ohi = [int(a.max()) for a in open_idx]
    while len(olo) < 3:        # 2D run: flatten the channel onto z=0
        olo.append(0)
        ohi.append(0)

    xs, ys, zs = _box_edges(olo, ohi)
    traces.append(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", name="Device channel", hoverinfo="skip",
        line=dict(color="#2dd4c4" if dark else "#0c8f96", width=3)))
    return traces


def _agent_traces(groups, theme="showcase", glow=True, marker_size=4):
    """Marker traces for one timepoint, in a fixed order (see ``AGENT_LABELS``).

    With ``glow`` each group emits a large translucent halo behind a crisp core,
    approximating the fluorescence look without a custom renderer.
    """
    import plotly.graph_objects as go

    style = THEMES[theme]
    traces = []
    for label in AGENT_LABELS:
        xs, ys, zs = groups.get(label, ([], [], []))
        colour, symbol = style[label]
        if glow:
            traces.append(go.Scatter3d(
                x=xs, y=ys, z=zs, mode="markers", name=label, legendgroup=label,
                showlegend=False, hoverinfo="skip",
                marker=dict(size=marker_size * 3.2, color=colour,
                            symbol="circle", opacity=0.10)))
        traces.append(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="markers", name=label, legendgroup=label,
            showlegend=True,
            marker=dict(size=marker_size, color=colour, symbol=symbol, opacity=0.95)))
    return traces


def showcase_figure(data, grid=None, steps=None, theme="showcase", glow=True,
                    show_device=True, config=None, marker_size=4,
                    frame_duration=400):
    """An animated, presentation-styled 3D view of a run's agents.

    Plays through every saved step with native plotly controls (play/pause and a
    frame slider), so the result is self-contained enough to export to a single
    HTML file. ``theme="plain"`` keeps the original light palette.
    """
    import plotly.graph_objects as go

    cells = data["cells"]
    grid = data["grids"][0] if grid is None else grid
    steps = list(data["steps"] if steps is None else steps)
    if not steps:
        raise ValueError("the simulation has no saved steps to animate")
    if config is None:
        config = metaspread.configs      # populated by load_simulation

    animated = _agent_traces(agent_groups(cells, steps[0], grid), theme, glow, marker_size)
    device = device_frame_traces(config, _spatial_shape(data), theme) if show_device else []

    fig = go.Figure(data=animated + device)
    # frames replace the agent traces by index; the device wireframe stays put
    fig.frames = [go.Frame(
        name=str(step), traces=list(range(len(animated))),
        data=_agent_traces(agent_groups(cells, step, grid), theme, glow, marker_size))
        for step in steps]

    dark = theme == "showcase"
    axis = dict(showbackground=False, showgrid=True,
                gridcolor="#1d2735" if dark else "#e6e6e6",
                zeroline=False, color="#8492a3" if dark else "#444444")
    play_args = {"frame": {"duration": frame_duration, "redraw": True},
                 "fromcurrent": True, "transition": {"duration": 0}}
    pause_args = {"frame": {"duration": 0, "redraw": False},
                  "mode": "immediate", "transition": {"duration": 0}}

    fig.update_layout(
        template="plotly_dark" if dark else "plotly_white",
        paper_bgcolor=_DARK_BG if dark else "white",
        title=dict(text=f"{os.path.basename(data['path'])} — grid {grid}"),
        height=720,
        legend=dict(itemsizing="constant"),
        scene=dict(xaxis=dict(title="x", **axis), yaxis=dict(title="y", **axis),
                   zaxis=dict(title="z", **axis), aspectmode="data",
                   bgcolor=_DARK_BG if dark else "white",
                   camera=dict(eye=dict(x=1.6, y=1.6, z=0.9))),
        updatemenus=[dict(
            type="buttons", direction="left", showactive=False,
            x=0.02, y=0, xanchor="left", yanchor="top", pad=dict(t=60),
            buttons=[dict(label="▶ Play", method="animate", args=[None, play_args]),
                     dict(label="❚❚ Pause", method="animate", args=[[None], pause_args])])],
        sliders=[dict(
            active=0, x=0.12, len=0.86, pad=dict(t=60),
            currentvalue=dict(prefix="step "),
            steps=[dict(label=str(step), method="animate",
                        args=[[str(step)], {"frame": {"duration": 0, "redraw": True},
                                            "mode": "immediate"}])
                   for step in steps])],
    )
    return fig


def export_showcase_html(simulation_path, out_path, grid=None, **kwargs):
    """Render a run's showcase view to a standalone HTML file.

    plotly.js is inlined, so the result opens offline and can be shared as a
    single file — the same property as a hand-built page, but with real data.
    """
    data = load_simulation(simulation_path)
    fig = showcase_figure(data, grid=grid, **kwargs)
    fig.write_html(out_path, include_plotlyjs="inline", full_html=True,
                   config={"displaylogo": False})
    return out_path
