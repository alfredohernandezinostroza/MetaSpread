"""Interactive 3D viewer for a MetaSpread simulation (marimo app).

Local use (needs the optional 'viz' extra: marimo + plotly)::

    pixi run -e viz marimo run metaspread/viewer3d_app.py -- --sim Simulations/<name>

All data loading and figure building live in ``metaspread.viewer3d``; this file
is only the reactive UI (step / grid / field / z sliders wired to those helpers).
"""
import marimo

app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import sys

    import marimo as mo

    from metaspread import viewer3d
    return mo, os, sys, viewer3d


@app.cell
def _(mo, sys):
    # simulation path from the CLI: `marimo run viewer3d_app.py -- --sim <path>`
    # (falls back to a trailing positional arg or the METASPREAD_SIM env var)
    import os as _os
    _args = mo.cli_args()
    sim_path = _args.get("sim") or _os.environ.get("METASPREAD_SIM", "")
    if not sim_path and len(sys.argv) > 1 and not sys.argv[-1].endswith(".py"):
        sim_path = sys.argv[-1]
    mo.stop(not sim_path,
            mo.md("**Pass a simulation folder:** `marimo run metaspread/viewer3d_app.py "
                  "-- --sim Simulations/<name>`"))
    mo.md(f"# MetaSpread 3D viewer\n`{sim_path}`")
    return (sim_path,)


@app.cell
def _(sim_path, viewer3d):
    data = viewer3d.load_simulation(sim_path)
    return (data,)


@app.cell
def _(data, mo):
    step = mo.ui.slider(start=0, stop=len(data["steps"]) - 1, value=len(data["steps"]) - 1,
                        label="step")
    grid = mo.ui.dropdown(options={str(g): g for g in data["grids"]},
                          value=str(data["grids"][0]), label="grid")
    field = mo.ui.dropdown(options=["(none)"] + data["fields"], value="(none)", label="field")
    z = mo.ui.slider(start=0, stop=max(data["gridsize_z"] - 1, 0), value=data["gridsize_z"] // 2,
                     label="z-slice")
    controls = mo.hstack([step, grid, field, z], justify="start")
    controls
    return field, grid, step, z


@app.cell
def _(data, grid, mo, step, viewer3d):
    current_step = data["steps"][step.value]
    groups = viewer3d.agent_groups(data["cells"], current_step, grid.value)
    fig = viewer3d.agent_scatter_figure(groups)
    fig.update_layout(title=f"Agents — step {current_step}, grid {grid.value}")
    mo.ui.plotly(fig)
    return (current_step,)


@app.cell
def _(current_step, field, grid, mo, sim_path, viewer3d, z):
    # optional field z-slice panel, shown only when a field is selected
    mo.stop(field.value == "(none)")
    slice2d = viewer3d.field_slice(sim_path, field.value, grid.value, current_step, z.value)
    mo.ui.plotly(viewer3d.field_heatmap_figure(
        slice2d, title=f"{field.value} — z={z.value}, step {current_step}, grid {grid.value}"))
    return


if __name__ == "__main__":
    app.run()
