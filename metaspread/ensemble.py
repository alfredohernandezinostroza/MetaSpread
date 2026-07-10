"""Run batches of simulations for replicates and parameter sweeps.

This is the entry point that later phases (calibration, sensitivity analysis,
ML-surrogate training) rely on: it runs many in-memory simulations while varying
the random seed and/or parameter values, and returns a single tidy DataFrame
tagged with ``run_id`` and the seed/overrides used.
"""
import pandas as pd

from metaspread import api


def run_ensemble(config=None, max_steps=None, data_collection_period=None,
                 seeds=None, param_overrides=None):
    """Run one simulation per (param_overrides x seeds) combination.

    Parameters
    ----------
    config : Config | str | None
        Base configuration (see :func:`metaspread.api.run`).
    max_steps, data_collection_period : int
        Passed through to each run.
    seeds : list[int] | None
        Random seeds to run. ``None`` means a single run with seed ``None``.
    param_overrides : list[dict] | None
        Each dict is applied via ``config.copy(**overrides)`` before running.
        ``None`` means a single run with the base config.

    Returns
    -------
    (combined_df, results)
        ``combined_df`` concatenates every run's ``agent_data`` with ``run_id``,
        ``seed`` and one column per overridden parameter. ``results`` is the list
        of :class:`SimulationResults` in run order.
    """
    base_config = api._coerce_config(config)
    seeds = [None] if seeds is None else list(seeds)
    param_overrides = [{}] if param_overrides is None else list(param_overrides)

    frames = []
    results = []
    run_id = 0
    for overrides in param_overrides:
        run_config = base_config.copy(**overrides) if overrides else base_config
        for seed in seeds:
            res = api.run(run_config, max_steps, data_collection_period,
                          seed=seed, save_path=None)
            df = res.agent_data.copy()
            df["run_id"] = run_id
            df["seed"] = seed
            for name, value in overrides.items():
                df[name] = str(value)
            frames.append(df)
            results.append(res)
            run_id += 1

    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return combined, results
