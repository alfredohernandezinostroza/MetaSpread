from numpy import number
import pandas as pd
import os
import ast

# Canonical ordered list of the simulation parameters stored in
# simulations_configs.csv. Keeping this in one place lets the Config object,
# the default-config generator and the validation logic stay in sync.
PARAM_NAMES = [
    "th", "tha", "xh", "xha", "dM", "dE", "phiM", "phiE", "dmmp", "theta",
    "Lambda", "gamma1", "gamma2", "vasculature_time", "doubling_time_M",
    "doubling_time_E", "single_cell_survival", "cluster_survival",
    "extravasation_probs", "dissagreggation_prob", "carrying_capacity",
    "normal_vessels_primary", "ruptured_vessels_primary",
    "secondary_sites_vessels", "n_center_points_for_tumor",
    "n_center_points_for_Vessels", "gridsize", "grids_number",
    "mesenchymal_proportion", "epithelial_proportion", "number_of_initial_cells",
]

# Extra keys that a saved simulation's configs.csv carries in addition to the
# core parameters above.
RUNTIME_NAMES = ["max_steps", "data_collection_period"]


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
        df = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
        return df

    @classmethod
    def from_csv(cls, path, validate=True):
        """Build a Config from a simulations_configs.csv (31 core params)."""
        df = cls._read_csv(path)
        values = dict(zip(df["Names"], df["Values"]))
        if validate:
            validate_configs(values)
        return cls(values)

    @classmethod
    def from_saved_simulation(cls, path):
        """Build a Config from a saved simulation's configs.csv.

        The saved file appends max_steps and data_collection_period after the
        core params; those are dropped here because they are provided per-run.
        """
        df = cls._read_csv(path)
        df = df[df["Names"].isin(PARAM_NAMES)]
        if len(df) != len(PARAM_NAMES):
            raise Exception(
                f"Expected {len(PARAM_NAMES)} configuration options, found {len(df)}!"
            )
        values = dict(zip(df["Names"], df["Values"]))
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
    df_configs = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
    if len(df_configs) < 33:
        raise Exception("Less than 33 configuration options! Are there some missing?")
    if len(df_configs) > 33:
        raise Exception("More than 33 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)
    return(list(df_configs["Names"]))

def load_simulation_configs_for_reloaded_simulation(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """
    df_configs = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
    # if the configs were loaded to continue from a previous
    # simulation, drop the last two rows
    # they shoud not be loaded as they change for each simulation according to user input
    # (max steps and step size)
    df_configs = df_configs[:-2]
    if len(df_configs) < 31:
        raise Exception("Less than 33 configuration options! Are there some missing?")
    if len(df_configs) > 31:
        raise Exception("More than 33 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)
    return(list(df_configs["Names"]))

def generate_default_configs():
    """Creates a default simulation_configs.csv file"""
    names = ["th","tha","xh","xha","dM","dE","phiM","phiE","dmmp","theta","Lambda","gamma1","gamma2","vasculature_time","doubling_time_M","doubling_time_E","single_cell_survival","cluster_survival","extravasation_probs","dissagreggation_prob","carrying_capacity","normal_vessels_primary","ruptured_vessels_primary","secondary_sites_vessels","n_center_points_for_tumor","n_center_points_for_Vessels","gridsize","grids_number","mesenchymal_proportion","epithelial_proportion","number_of_initial_cells"]
    values = [0.001,0.001,0.005,0.005,1e-4,5e-5,0.0005,0.0005,0.001,0.195,0.1,1,1,180,2000,3000,5e-04,0.025,[0.75, 0.25],0.5,4,8,2,[10, 10],97,200,201,3,0.6,0.4,388]
    default_configs = pd.DataFrame({"Names": names, "Values": values})
    default_configs.to_csv("simulations_configs.csv", index=False)

def check_if_configs_are_present():
    """Checks if a configs file is present in the current directory"""
    if not os.path.isfile("simulations_configs.csv"):
        print("Configs file`for new simulation \"simulations_configs.csv\" not found! Creating a configs file with default values.")
        generate_default_configs()
