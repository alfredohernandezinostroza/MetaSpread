from numpy import number
import pandas as pd
import os
import ast


def _safe_literal(value):
    """Parse a configs.csv cell value.

    Values are stored as Python literals (numbers, bools, lists) and read back
    with ast.literal_eval. String-valued parameters (e.g. device_mask_path) may
    be written as a bare string or left empty; a value that is not a valid literal
    (including an empty cell) falls back to the raw string. Every pre-existing
    parameter is already a valid literal, so this changes none of their parsing.
    """
    if isinstance(value, str):
        if value.strip() == "":
            return ""
        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return value
    return value

# Canonical ordered list of the simulation parameters stored in
# simulations_configs.csv. Keeping this in one place lets the Config object,
# the default-config generator and the validation logic stay in sync.
# The original (Phase 0) core parameters, followed by the Phase 1 additions.
# Phase 1 parameters are all inert at their defaults so a default-config run is
# byte-identical to the pre-Phase-1 output.
_CORE_PARAM_NAMES = [
    "th", "tha", "xh", "xha", "dM", "dE", "phiM", "phiE", "dmmp", "theta",
    "Lambda", "gamma1", "gamma2", "vasculature_time", "doubling_time_M",
    "doubling_time_E", "single_cell_survival", "cluster_survival",
    "extravasation_probs", "dissagreggation_prob", "carrying_capacity",
    "normal_vessels_primary", "ruptured_vessels_primary",
    "secondary_sites_vessels", "n_center_points_for_tumor",
    "n_center_points_for_Vessels", "gridsize", "grids_number",
    "mesenchymal_proportion", "epithelial_proportion", "number_of_initial_cells",
]

# Phase 1 — EMT/MET plasticity, oxygen field, immune agents (all off by default).
_PHASE1_PARAM_NAMES = [
    # EMT/MET plasticity
    "emt_prob", "met_prob", "enable_hypoxia_emt", "hypoxia_threshold",
    # oxygen / nutrient field
    "enable_oxygen", "d_oxygen", "oxygen_supply", "oxygen_consumption",
    "oxygen_initial", "oxygen_max",
    # immune agents
    "enable_immune", "n_immune_cells", "immune_kill_prob", "immune_diff_coeff",
]

# Phase 3D — optional 3D spatial domain. space_dimensions defaults to 2 (the
# existing 2D path, byte-identical); gridsize_z is only used when it is 3.
_PHASE3D_PARAM_NAMES = ["space_dimensions", "gridsize_z"]

# Phase 2 — organ-on-chip (all off by default, so a default run is byte-identical).
_PHASE2_PARAM_NAMES = [
    # flow / advection of the diffusible fields
    "enable_flow", "flow_velocity",
    # shear-dependent survival of circulating clusters
    "enable_shear", "shear_stress", "shear_death_coeff",
    # device geometry: impassable walls (parametric channel or a mask file)
    "enable_device_geometry", "channel_axis", "channel_margin", "device_mask_path",
]

# Canonical ordered list of every simulation parameter stored in
# simulations_configs.csv. Keeping this in one place lets the Config object, the
# default-config generator and the validation logic stay in sync.
PARAM_NAMES = (
    _CORE_PARAM_NAMES + _PHASE1_PARAM_NAMES + _PHASE3D_PARAM_NAMES
    + _PHASE2_PARAM_NAMES
)

# Extra keys that a saved simulation's configs.csv carries in addition to the
# core parameters above.
RUNTIME_NAMES = ["max_steps", "data_collection_period"]

# Single source of truth for default parameter values. Used to generate the
# default simulations_configs.csv and to backfill parameters that are missing
# from an older saved simulation's configs.csv.
DEFAULTS = {
    "th": 0.001, "tha": 0.001, "xh": 0.005, "xha": 0.005,
    "dM": 1e-4, "dE": 5e-5, "phiM": 0.0005, "phiE": 0.0005, "dmmp": 0.001,
    "theta": 0.195, "Lambda": 0.1, "gamma1": 1, "gamma2": 1,
    "vasculature_time": 180, "doubling_time_M": 2000, "doubling_time_E": 3000,
    "single_cell_survival": 5e-04, "cluster_survival": 0.025,
    "extravasation_probs": [0.75, 0.25], "dissagreggation_prob": 0.5,
    "carrying_capacity": 4, "normal_vessels_primary": 8,
    "ruptured_vessels_primary": 2, "secondary_sites_vessels": [10, 10],
    "n_center_points_for_tumor": 97, "n_center_points_for_Vessels": 200,
    "gridsize": 201, "grids_number": 3,
    "mesenchymal_proportion": 0.6, "epithelial_proportion": 0.4,
    "number_of_initial_cells": 388,
    # --- Phase 1 (inert defaults) ---
    "emt_prob": 0.0, "met_prob": 0.0,
    "enable_hypoxia_emt": False, "hypoxia_threshold": 0.1,
    "enable_oxygen": False, "d_oxygen": 0.001, "oxygen_supply": 0.1,
    "oxygen_consumption": 0.01, "oxygen_initial": 1.0, "oxygen_max": 1.0,
    "enable_immune": False, "n_immune_cells": 0,
    "immune_kill_prob": 0.1, "immune_diff_coeff": 1e-4,
    # --- Phase 3D (2D by default) ---
    "space_dimensions": 2, "gridsize_z": 41,
    # --- Phase 2 (organ-on-chip; inert defaults) ---
    # flow_velocity is a per-axis advection velocity; only the first
    # space_dimensions entries are used, so one default serves 2D and 3D.
    "enable_flow": False, "flow_velocity": [0.0, 0.0, 0.0],
    "enable_shear": False, "shear_stress": 0.0, "shear_death_coeff": 0.0,
    "enable_device_geometry": False, "channel_axis": 0, "channel_margin": 0,
    "device_mask_path": "",
}


def _backfill_defaults(values):
    """Fill any missing PARAM_NAMES in `values` from DEFAULTS (in place).

    Lets configs.csv files written before a parameter existed (e.g. pre-Phase-1
    simulations) still load, with the new parameters taking their inert default.
    """
    import warnings
    missing = [name for name in PARAM_NAMES if name not in values]
    if missing:
        warnings.warn(
            f"Config is missing parameters {missing}; using defaults for them."
        )
        for name in missing:
            values[name] = DEFAULTS[name]
    return values


def validate_configs(d):
    """Validate a dict of config values. Raises ValueError with all problems.

    This is the same set of checks that historically lived inside
    init_simulation_configs, extracted so both the legacy functions and the
    Config object can share it.
    """
    error_string = ""
    if sum(d["extravasation_probs"]) != 1:
        error_string += "Extravasation probabilities must sum 1!\n"
    if len(d["extravasation_probs"]) != d["grids_number"] - 1:
        error_string += "There must be as many Extravasation probabilities as the value of (grids_number - 1)!\n"
    if len(d["secondary_sites_vessels"]) != d["grids_number"] - 1:
        error_string += "There must be as many secondary site vessels as the value of grids_number - 1!\n"
    if d["mesenchymal_proportion"] + d["epithelial_proportion"] != 1:
        error_string += "Mesenchymal_proportion + epithelial_proportion must be 1!\n"
    if d["n_center_points_for_tumor"] <= 0:
        error_string += "n_center_points_for_tumor must be greater than 0!\n"
    if d["n_center_points_for_tumor"] > d["gridsize"]:
        error_string += "n_center_points_for_tumor must be less than or equal to gridsize!\n"
    if d["number_of_initial_cells"] <= 0:
        error_string += "number_of_initial_cells must be greater than 0!\n"
    if d["number_of_initial_cells"] > d["n_center_points_for_tumor"] * d["carrying_capacity"]:
        error_string += (
            f"number_of_initial_cells ({d['number_of_initial_cells']}) must be less than or equal to "
            f"n_center_points_for_tumor * carrying_capacity ({d['n_center_points_for_tumor'] * d['carrying_capacity']})!\n"
        )

    # --- Phase 1 checks (only bite when the relevant feature is enabled) ---
    for prob_name in ("emt_prob", "met_prob", "immune_kill_prob"):
        prob = d.get(prob_name, 0.0)
        if not (0.0 <= prob <= 1.0):
            error_string += f"{prob_name} must be between 0 and 1!\n"
    if d.get("enable_hypoxia_emt", False) and not d.get("enable_oxygen", False):
        error_string += "enable_hypoxia_emt requires enable_oxygen to be True!\n"
    if d.get("n_immune_cells", 0) < 0:
        error_string += "n_immune_cells must be >= 0!\n"
    if d.get("space_dimensions", 2) not in (2, 3):
        error_string += "space_dimensions must be 2 or 3!\n"
    if d.get("space_dimensions", 2) == 3 and d.get("gridsize_z", 1) <= 0:
        error_string += "gridsize_z must be greater than 0 when space_dimensions == 3!\n"

    # --- Phase 2 checks (only bite when the relevant feature is enabled) ---
    if d.get("enable_flow", False):
        ndim = d.get("space_dimensions", 2)
        fv = d.get("flow_velocity", [])
        if not isinstance(fv, (list, tuple)) or len(fv) < ndim:
            error_string += (
                f"flow_velocity must be a list of at least space_dimensions ({ndim}) "
                "numbers when enable_flow is True!\n"
            )
        elif not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in fv):
            error_string += "flow_velocity entries must be numbers!\n"
    if d.get("enable_shear", False):
        # both non-negative keeps the survival factor exp(-coeff*stress) in (0, 1]
        if d.get("shear_stress", 0.0) < 0:
            error_string += "shear_stress must be >= 0 when enable_shear is True!\n"
        if d.get("shear_death_coeff", 0.0) < 0:
            error_string += "shear_death_coeff must be >= 0 when enable_shear is True!\n"
    if d.get("enable_device_geometry", False) and not d.get("device_mask_path", ""):
        ndim = d.get("space_dimensions", 2)
        axis = d.get("channel_axis", 0)
        if not (isinstance(axis, int) and 0 <= axis < ndim):
            error_string += f"channel_axis must be an integer in [0, space_dimensions) ({ndim})!\n"
        if d.get("channel_margin", 0) < 0:
            error_string += "channel_margin must be >= 0!\n"

    if error_string != "":
        raise ValueError(error_string)


class Config:
    """Holds a full set of simulation parameters as attributes.

    Replaces the previous pattern of injecting parameters as module globals,
    which limited the process to a single parameter set. A Config instance can
    be passed to CancerModel, copied with overrides for parameter sweeps, and
    saved/loaded to CSV. For backward compatibility with modules and tests that
    still read ``metaspread.configs.<name>`` (postprocessing, legacy tests),
    ``publish_to_module`` mirrors the values onto this module's namespace.
    """

    def __init__(self, values):
        # values: dict of name -> value (may include RUNTIME_NAMES)
        self._names = [n for n in PARAM_NAMES if n in values]
        for name, value in values.items():
            setattr(self, name, value)

    @property
    def param_names(self):
        """Names of the core simulation parameters this Config carries."""
        return list(self._names)

    @classmethod
    def _read_csv(cls, path):
        df = pd.read_csv(path, header=0, converters={"Values": _safe_literal})
        return df

    @classmethod
    def from_csv(cls, path, validate=True):
        """Build a Config from a simulations_configs.csv.

        Parameters missing from the file (e.g. Phase 1 params in an older config)
        are backfilled from DEFAULTS.
        """
        df = cls._read_csv(path)
        values = dict(zip(df["Names"], df["Values"]))
        _backfill_defaults(values)
        if validate:
            validate_configs(values)
        return cls(values)

    @classmethod
    def from_saved_simulation(cls, path):
        """Build a Config from a saved simulation's configs.csv.

        The saved file appends max_steps and data_collection_period after the
        core params; those are dropped here because they are provided per-run.
        Parameters absent from an older saved simulation are backfilled from
        DEFAULTS so pre-Phase-1 simulations still load.
        """
        df = cls._read_csv(path)
        df = df[df["Names"].isin(PARAM_NAMES)]
        values = dict(zip(df["Names"], df["Values"]))
        _backfill_defaults(values)
        return cls(values)

    def copy(self, **overrides):
        """Return a new Config with the given parameter values overridden."""
        values = {name: getattr(self, name) for name in self._names}
        values.update(overrides)
        return Config(values)

    def as_dict(self):
        return {name: getattr(self, name) for name in self._names}

    def to_csv(self, path, extra=None):
        """Write this Config (plus optional extra name->value) to CSV."""
        names = list(self._names)
        values = [getattr(self, name) for name in names]
        if extra:
            for name, value in extra.items():
                names.append(name)
                values.append(value)
        df_vars = pd.DataFrame({"Names": names, "Values": values}).set_index("Names")
        df_vars.to_csv(path)

    def publish_to_module(self):
        """Mirror values onto metaspread.configs.<name> for backward compat."""
        globals().update(self.as_dict())


def init_simulation_configs(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """
    config = Config.from_csv(path, validate=True)
    config.publish_to_module()
    return config.param_names

def load_simulation_configs_for_data_generation(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """
    df_configs = pd.read_csv(path, header=0, converters={"Values": _safe_literal})
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    for rt in RUNTIME_NAMES:
        if rt not in dict_configs:
            raise Exception(f"Missing runtime parameter '{rt}' in {path}!")
    _backfill_defaults(dict_configs)
    globals().update(dict_configs)
    return(list(dict_configs.keys()))

def load_simulation_configs_for_reloaded_simulation(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """
    df_configs = pd.read_csv(path, header=0, converters={"Values": _safe_literal})
    # drop the runtime rows (max_steps, data_collection_period): they change for
    # each simulation according to user input
    df_configs = df_configs[df_configs["Names"].isin(PARAM_NAMES)]
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    _backfill_defaults(dict_configs)
    globals().update(dict_configs)
    return(list(dict_configs.keys()))

def generate_default_configs():
    """Creates a default simulation_configs.csv file"""
    names = list(PARAM_NAMES)
    values = [DEFAULTS[name] for name in names]
    default_configs = pd.DataFrame({"Names": names, "Values": values})
    default_configs.to_csv("simulations_configs.csv", index=False)

def check_if_configs_are_present():
    """Checks if a configs file is present in the current directory"""
    if not os.path.isfile("simulations_configs.csv"):
        print("Configs file`for new simulation \"simulations_configs.csv\" not found! Creating a configs file with default values.")
        generate_default_configs()
