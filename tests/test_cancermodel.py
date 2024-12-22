import pytest
from metaspread import cancermodel
from metaspread import configs

def test_imports():
    assert cancermodel.mesa is not None
    assert cancermodel.plt is not None
    assert cancermodel.np is not None
    assert cancermodel.pd is not None
    assert cancermodel.os is not None
    assert cancermodel.json is not None
    assert cancermodel.ast is not None
    assert cancermodel.CancerCell is not None
    assert cancermodel.Vessel is not None
    assert cancermodel.find_quasi_circle is not None

# def test_generate_cancer_model(mocker):
#     mocker.patch('metaspread.cancermodel.generate_cancer_model')
#     cancermodel.generate_cancer_model()
#     assert cancermodel.generate_cancer_model.called

from metaspread.cancermodel import CancerModel
from metaspread.cancercell import CancerCell
from metaspread.vessel import Vessel
import numpy as np
import pandas as pd
import pytest
import ast

#todo: model is not callable (duh! I think I cannot call a private variable (is it though?))
#todo: use tmp_path_facorty to create the model once, and use it for the rest of the tests

def test_cancermodel(tmp_path) -> None:
    temp_simulation_folder = tmp_path / "Simulations" 
    temp_simulation_folder.mkdir()
    model = CancerModel(
        number_of_initial_cells=30,
        width=201,
        height=201,
        grids_number=3,
        max_steps=1000,
        data_collection_period=200000,
        new_simulation_folder=temp_simulation_folder)
    
    assert model.number_of_initial_cells==30
    assert model.width==configs.gridsize
    assert model.height==configs.gridsize
    assert model.grids_number==configs.grids_number
    assert model.max_steps==1000
    assert model.data_collection_period==200000
    assert model.new_simulation_folder==temp_simulation_folder

    #test proliferation
    cancermodel.carrying_capacity = 200
    # current_cell_count = list(map(type, model.schedule.agents)).count(CancerCell)
    current_cell_count = cancermodel.count_total_cells(model)
    model.proliferate("mesenchymal")
    model.proliferate("epithelial")
    assert 2*current_cell_count == list(map(type, model.schedule.agents)).count(CancerCell)
    current_cell_count = 2*current_cell_count

    #test calculate_environment
    grids_number = model.grids_number
    for i in range(grids_number):
        assert (model.ecm[i]  == 1.0).all() == True
        assert (model.mmp2[i] == 0.0).all() == True
    model.calculate_environment(model.mmp2, model.ecm)
    assert (model.ecm[0]  == 1.0).all() == False
    assert (model.mmp2[0] == 0.0).all() == False
    for i in range(1,grids_number):
        assert (model.ecm[i]  == 1.0).all() == True 
        assert (model.mmp2[i] == 0.0).all() == True
    
    #test count_vasculature_cells
    assert cancermodel.count_vasculature_cells(model) == 0
    model.vasculature = {1: [(100,100)]}
    assert cancermodel.count_vasculature_cells(model) == 200

    #test cell disaggregation
    for i in range(100):
        for j in range(100):
            model.vasculature = {1: [(i,j)]}
            model.disaggregate_clusters(1)
            assert sum(x+y for x,y in model.vasculature[1]) == i + j

    #test cell travel
    model.vasculature = {1: [(100,100)]}
    cancermodel.single_cell_survival = 1
    cancermodel.cluster_survival = 1
    model.step()
    model.step()
    assert list(map(type, model.schedule.agents)).count(CancerCell) == current_cell_count + 200
    current_cell_count = current_cell_count + 200
    assert current_cell_count == list(map(type, model.schedule.agents)).count(CancerCell)
    

    
