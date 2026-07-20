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
