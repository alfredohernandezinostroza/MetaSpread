"""Tests for Phase 5: visualization/postprocessing of the new simulation outputs.

Postprocessing is read-only over saved artifacts, so it cannot affect the
byte-identity sim gate. The contract these tests pin down is:
  * an oxygen-enabled run produces Oxygen graphs (and videos-ready PNGs),
  * a run without oxygen produces no oxygen artifacts at all (off by default),
  * the amount_of_pictures==0 ("all pictures") path renders every field/grid,
    not just the first (regression guard for the range_of_pictures fix).
"""
import matplotlib
matplotlib.use("Agg")  # headless: never try to open a window during tests

from metaspread import Config, run
from metaspread import datagenerator, graphgenerator


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


def _run_and_postprocess(cfg, tmp_path, monkeypatch, max_steps=4, period=2):
    run(cfg, max_steps, period, seed=1, save_path=tmp_path)
    name = _sim_name(cfg, max_steps, period)
    monkeypatch.chdir(tmp_path)          # generators resolve "Simulations/<name>" from CWD
    datagenerator.generate_data(name)
    graphgenerator.generate_graphs(name, 0)   # 0 == "all pictures"
    return tmp_path / "Simulations" / name


def test_oxygen_run_produces_oxygen_graphs(tmp_path, monkeypatch):
    cfg = _small_config(enable_oxygen=True)
    sim_dir = _run_and_postprocess(cfg, tmp_path, monkeypatch)

    oxygen_images = sim_dir / "Graphical analysis" / "Oxygen dynamics"
    assert oxygen_images.is_dir()
    pngs = list(oxygen_images.glob("Oxygen-grid*-step*.png"))
    assert pngs, "no oxygen images were rendered"
    # one image per (grid, saved step): 2 grids x steps {2,4}
    assert len(pngs) == cfg.grids_number * 2


def test_no_oxygen_run_produces_no_oxygen_artifacts(tmp_path, monkeypatch):
    cfg = _small_config()  # oxygen off (default)
    sim_dir = _run_and_postprocess(cfg, tmp_path, monkeypatch)

    assert not (sim_dir / "Oxygen").exists()                       # no raw field
    assert not (sim_dir / "Graphical analysis" / "Oxygen dynamics").exists()  # no images


def test_immune_run_marks_immune_cells_in_tumor_coords(tmp_path, monkeypatch):
    import pandas as pd
    # kill_prob 0 so the immune agents persist and show up at every saved step
    cfg = _small_config(enable_immune=True, n_immune_cells=15, immune_kill_prob=0.0)
    sim_dir = _run_and_postprocess(cfg, tmp_path, monkeypatch)

    tumor_data = sim_dir / "Data analysis" / "Tumor dynamics"
    coords_files = sorted(tumor_data.glob("Cells-grid1-step*Tumor size at*.csv"))
    assert coords_files, "no tumor coords files were written"
    coords = pd.read_csv(coords_files[-1], index_col=0)
    assert len(coords) == 10                 # immune X/Y appended as rows 8-9
    assert coords.iloc[8].notna().any()      # grid 1 actually carries immune positions
    # and the scatter (which reads those rows) rendered without error
    assert list((sim_dir / "Graphical analysis" / "Tumor dynamics").glob("*grid1-step*.png"))


def test_3d_data_generalizes_histogram_and_centroid(tmp_path, monkeypatch):
    # 3D fields are saved as .npy, which used to make generate_data early-return.
    # It must now run and produce n-D-correct analytics (histogram over the full
    # gridsize^2*gridsize_z lattice; a centroid with a z coordinate).
    import pandas as pd
    from metaspread import datagenerator
    cfg = _small_config(space_dimensions=3, gridsize=11, gridsize_z=5,
                        number_of_initial_cells=8, n_center_points_for_tumor=8)
    run(cfg, 4, 2, seed=1, save_path=tmp_path)
    name = _sim_name(cfg, 4, 2)
    monkeypatch.chdir(tmp_path)
    datagenerator.generate_data(name)

    tumor = tmp_path / "Simulations" / name / "Data analysis" / "Tumor dynamics"
    hist_files = sorted(tumor.glob("*Histogram*.csv"))
    assert hist_files, "generate_data early-returned for the 3D run"
    h = pd.read_csv(hist_files[-1], index_col=0)
    assert int(h["Frequency"].sum()) == cfg.gridsize * cfg.gridsize * cfg.gridsize_z
    rad = pd.read_csv(tumor / "Tumor radius and diameter history in grid 1.csv", index_col=0)
    assert "Centroid z" in rad.columns
    assert rad["Radius"].notna().any()


def test_3d_run_renders_field_montages(tmp_path, monkeypatch):
    # In 3D each field is a .npy volume; generate_graphs must render it as a
    # multi-panel z-slice montage instead of crashing on the missing CSVs.
    import matplotlib.pyplot as plt
    cfg = _small_config(space_dimensions=3, gridsize=11, gridsize_z=5,
                        number_of_initial_cells=8, n_center_points_for_tumor=8,
                        enable_oxygen=True)
    sim_dir = _run_and_postprocess(cfg, tmp_path, monkeypatch)

    ga = sim_dir / "Graphical analysis"
    for sub in ("Ecm dynamics", "Mmp2 dynamics", "Oxygen dynamics"):
        assert list((ga / sub).glob("*grid1-step*.png")), f"{sub} produced no 3D montage"
    # a z-slice montage is multiple panels wide, unlike the ~600px single 2D heatmap
    montage = sorted((ga / "Ecm dynamics").glob("*.png"))[0]
    assert plt.imread(montage).shape[1] >= 800


def test_all_pictures_renders_every_field_and_grid(tmp_path, monkeypatch):
    # Regression for the range_of_pictures exhaustion bug: with amount==0 the
    # generator must render Mmp2, Ecm and Tumor for BOTH grids, not just grid 1.
    cfg = _small_config()
    sim_dir = _run_and_postprocess(cfg, tmp_path, monkeypatch)
    graphical = sim_dir / "Graphical analysis"
    for subdir in ("Mmp2 dynamics", "Ecm dynamics", "Tumor dynamics"):
        for grid_id in range(1, cfg.grids_number + 1):
            hits = list((graphical / subdir).glob(f"*grid{grid_id}-step*.png"))
            assert hits, f"{subdir} produced nothing for grid {grid_id}"
