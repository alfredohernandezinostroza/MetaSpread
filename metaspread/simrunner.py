import metaspread.configs
import pandas as pd
import shutil
import mesa
import ast
import os
from pathlib import Path

# To run this code you must be in the parent folder of the program


def run_simulation(max_steps, data_collection_period, save_path=Path("."), loaded_simulation_path="", config=None, seed=None):

    # load configs file from a previous simulation or loads the general configs file
    print(loaded_simulation_path)
    loaded_simulation_path = str(loaded_simulation_path).strip('\"') if loaded_simulation_path else ""
    if config is None:
        if loaded_simulation_path != "":
            configs_path = os.path.join(loaded_simulation_path, "configs.csv")
            config = metaspread.configs.Config.from_saved_simulation(configs_path)
        else:
            config = metaspread.configs.Config.from_csv("simulations_configs.csv")
    # mirror onto the configs module for backward-compat consumers
    config.publish_to_module()

    # Parameters for this simulation
    number_of_initial_cells = config.number_of_initial_cells  # Number of cancer cells
    gridsize     = config.gridsize
    grids_number = config.grids_number
    width        = gridsize
    height       = gridsize

    # Name of the directories
    simulations_dir = save_path / "Simulations"
    os.makedirs(simulations_dir, exist_ok=True)
    if loaded_simulation_path != "":
        cells_path = os.path.join(loaded_simulation_path, "CellsData.csv")
        df = pd.read_csv(cells_path)
        loaded_max_step = max(df["Step"])
        new_simulation_folder = os.path.normpath(loaded_simulation_path)
        new_simulation_folder = os.path.basename(new_simulation_folder)
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
    else:
        df = pd.DataFrame()
        loaded_max_step = 0
        new_simulation_folder = f"Sim-max_steps-{max_steps}-collection_period-{data_collection_period}-cells-{number_of_initial_cells}-grids_number-{grids_number}"

        # Creates the path for the new simulation
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
        pathMmp2 = os.path.join(new_simulation_path, "Mmp2")
        pathEcm = os.path.join(new_simulation_path, "Ecm")
        pathVasculature = os.path.join(new_simulation_path, "Vasculature")
        pathTimeOfPopulation = os.path.join(new_simulation_path, "Time when grids were populated")

        # Create folder for all cells analysis, for Mmp2 matrices and Ecm matrices
        if not os.path.exists(new_simulation_path):
            print(f'\t Folder for this simulation: {new_simulation_path}')
            print(f'\t Saving agents data at: {new_simulation_path}')
            print(f'\t Saving Mmp2 data at: {pathMmp2}')
            print(f'\t Saving Ecm data at: {pathEcm}')
            print(f'\t Saving Vasculature data at: {pathVasculature}')

            os.makedirs(new_simulation_path)
            os.makedirs(pathMmp2)
            os.makedirs(pathEcm)
            os.makedirs(pathVasculature)
            os.makedirs(pathTimeOfPopulation)
        # If there is already a simulation you skip it
        else:
            return print("This simulation already exists!")

    # Run the simulation and saves the data
    configs_save_path = os.path.join(new_simulation_path, 'configs.csv')
    print(f"\t Saving all the simulations parameters at: {configs_save_path}")
    config.to_csv(configs_save_path, extra={"max_steps": max_steps, "data_collection_period": data_collection_period})
    model = metaspread.CancerModel(
        number_of_initial_cells,
        width,
        height,
        grids_number,
        max_steps,
        data_collection_period,
        new_simulation_path,
        loaded_simulation_path,
        seed=seed,
        config=config)
    for i in range(max_steps):
        model.step()
    print(f'Finished the simulation at time step {model.schedule.time}!')
    return model
