from metaspread.cancermodel import CancerModel
from metaspread.cancercell import CancerCell
from metaspread.vessel import Vessel
import numpy as np
import pandas as pd
import pytest
import ast


def test_carrying_capacity_guard_with_vessel(tmp_path) -> None:
    """Regression test for the CancerCell.move carrying-capacity guard (6d43512).

    The guard must count only cancer cells at the destination. The old bug used a
    leftover loop variable (`agent`) instead of the comprehension variable, which
    miscounts when a vessel shares the destination cell. We fill a destination
    that already holds a vessel, one mover at a time, on a secondary grid (so the
    intravasation branch does not trigger), and assert the cell occupancy settles
    at exactly carrying_capacity.
    """
    folder = tmp_path / "cap_guard"
    folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0, width=201, height=201, grids_number=2,
        max_steps=10, data_collection_period=10, new_simulation_folder=folder,
        fixed_p_left=0, fixed_p_right=1, fixed_p_top=0, fixed_p_bottom=0,
    )
    grid_id = 2               # secondary grid -> no intravasation branch
    grid = model.grids[grid_id - 1]
    cap = model.config.carrying_capacity
    src = (100, 100)
    dest = (101, 100)

    # a vessel occupies the destination cell (does not count toward capacity)
    vessel = Vessel(model.current_agent_id, model, False, grid, grid_id)
    model.current_agent_id += 1
    grid.place_agent(vessel, dest)

    # more movers than capacity, all forced to step right into `dest`
    movers = []
    for _ in range(cap + 2):
        m = CancerCell(model.current_agent_id, model, grid, grid_id,
                       "mesenchymal", model.ecm[grid_id - 1], model.mmp2[grid_id - 1])
        model.current_agent_id += 1
        grid.place_agent(m, src)
        model.schedule.add(m)
        movers.append(m)

    for m in movers:
        if m.pos == src:      # only those still at the source can move
            m.move()

    cells_at_dest = len([a for a in grid.get_cell_list_contents([dest])
                         if a.agent_type == "cell"])
    assert cells_at_dest == cap

#todo: model is not callable (duh! I think I cannot call a private variable (is it though?))
#todo: use tmp_path_facorty to create the model once, and use it for the rest of the tests

def test_cancercell(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=30,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder
        )
    
    assert model.data_collection_period == 10
    assert model.number_of_initial_cells==30
    assert model.width==201
    assert model.height==201
    assert model.grids_number==2
    assert model.max_steps==1000
    assert model.data_collection_period==10
    assert model.new_simulation_folder==temp_simulation_folder

    ccell_id = None
    grid_id = 1
    #add to the following line the keyword arguments for the cancer cell
    ccell = CancerCell(
        unique_id=ccell_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )

    assert ccell.unique_id          == ccell_id
    assert ccell.model              == model
    assert ccell.grid               == model.grids[grid_id-1]
    assert ccell.grid_id            == grid_id
    assert ccell.phenotype          == "mesenchymal"
    assert np.array_equal(ccell.ecm, model.ecm[grid_id-1])
    assert np.array_equal(ccell.mmp2, model.mmp2[grid_id-1])
def test_cancercell_movement_left(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=1,
        fixed_p_right=0,
        fixed_p_top=0,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0]-1,pos[1]) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_right(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=1,
        fixed_p_top=0,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0]+1,pos[1]) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_down(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=0,
        fixed_p_top=0,
        fixed_p_bottom=1
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0],pos[1]-1) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_up(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=0,
        fixed_p_top=1,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        future_positions = [(pos[0],pos[1]+1) for pos in current_positions]
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == future_positions[i]
                i = i + 1

def test_cancercell_movement_stay(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "test_simulation_movement"
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=0,
        width=201,
        height=201,
        grids_number=2,
        max_steps=1000,
        data_collection_period=10,
        new_simulation_folder=temp_simulation_folder,
        fixed_p_left=0,
        fixed_p_right=0,
        fixed_p_top=0,
        fixed_p_bottom=0
        )

    grid_id = 1
    ccell = CancerCell(
        unique_id=model.current_agent_id,
        model=model, 
        grid=model.grids[grid_id-1], 
        grid_id=grid_id, 
        phenotype="mesenchymal", 
        ecm=model.ecm[grid_id-1], 
        mmp2=model.mmp2[grid_id-1]
    )
    x = 100
    y = 100
    model.current_agent_id += 1
    model.grids[grid_id-1].place_agent(ccell, (x,y)) 
    model.cancer_cells_counter[grid_id-1] += 1
    model.schedule.add(ccell)
    
    for j in range(1,10):
        current_positions = []
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                current_positions.append(agent.pos)
        model.step()
        i = 0
        for agent in model.schedule.agents:
            if agent.agent_type == "cell":
                assert agent.pos == current_positions[i]
                i = i + 1