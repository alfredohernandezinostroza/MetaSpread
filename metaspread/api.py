"""Programmatic entry point for running MetaSpread simulations.

The command-line interface (``python -m metaspread run ...``) writes every
artifact to disk. This module adds an in-process API that returns results as
pandas DataFrames, which is what parameter sweeps, calibration and ML-surrogate
training need. It is the foundation the ensemble runner (and later phases) build
on.

Typical use::

    import metaspread
    results = metaspread.run(max_steps=300, data_collection_period=30, seed=42)
    df = results.agent_data          # in-memory, nothing written to disk

    # or with an explicit / modified configuration
    cfg = metaspread.Config.from_csv("simulations_configs.csv")
    results = metaspread.run(cfg.copy(carrying_capacity=6), 300, 30, seed=1)
"""
import os
from pathlib import Path

import metaspread.configs
from metaspread.cancermodel import CancerModel


def _coerce_config(config):
    """Accept a Config, a path to a configs CSV, or None (load the default)."""
    if config is None:
        return metaspread.configs.Config.from_csv("simulations_configs.csv")
    if isinstance(config, metaspread.configs.Config):
        return config
    if isinstance(config, (str, os.PathLike)):
        return metaspread.configs.Config.from_csv(config)
    raise TypeError(
        f"config must be a Config, a path, or None; got {type(config).__name__}"
    )


class SimulationResults:
    """Lightweight accessor over a finished CancerModel's collected data."""

    def __init__(self, model):
        self.model = model
        self.config = model.config

    @property
    def agent_data(self):
        """Per-agent, per-step dataframe (Step, AgentID, Position, ...)."""
        df = self.model.datacollector.get_agent_vars_dataframe()
        return df.reset_index(level=["Step", "AgentID"])

    @property
    def model_data(self):
        """Per-step model-level dataframe (e.g. Total cells)."""
        return self.model.datacollector.get_model_vars_dataframe()


def run(config=None, max_steps=None, data_collection_period=None, seed=None,
        save_path=None, loaded_simulation_path=""):
    """Run a single simulation and return a :class:`SimulationResults`.

    Parameters
    ----------
    config : Config | str | os.PathLike | None
        Configuration to use. ``None`` loads ``simulations_configs.csv``.
    max_steps, data_collection_period : int
        Simulation length and how often data is collected.
    seed : int | None
        Random seed for reproducibility.
    save_path : str | os.PathLike | None
        If ``None`` (default), the simulation runs fully in memory and writes
        nothing to disk. Otherwise a ``Simulations/`` tree is created under this
        path exactly like the CLI.
    """
    if max_steps is None or data_collection_period is None:
        raise ValueError("max_steps and data_collection_period are required")

    config = _coerce_config(config)

    if save_path is None:
        model = CancerModel(
            config.number_of_initial_cells,
            config.gridsize,
            config.gridsize,
            config.grids_number,
            max_steps,
            data_collection_period,
            new_simulation_folder="",
            loaded_simulation_path=loaded_simulation_path,
            seed=seed,
            config=config,
            save_to_disk=False,
        )
        for _ in range(max_steps):
            model.step()
        print(f"\nFinished in-memory simulation at time step {model.schedule.time}!")
        return SimulationResults(model)

    # disk-backed path: reuse the CLI runner so folder layout stays identical
    from metaspread import simrunner
    model = simrunner.run_simulation(
        max_steps,
        data_collection_period,
        save_path=Path(save_path),
        loaded_simulation_path=loaded_simulation_path,
        config=config,
        seed=seed,
    )
    return SimulationResults(model)
