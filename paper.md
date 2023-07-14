---
title: 'MetaSpread: A cancer growth and metastatic spread simulation package in Python'
tags:
  - Python
  - agent based modelling
  - dynamics
  - partial differential equations
  - cancer
  - metastasis
authors:
  - name: Alfredo Hernández-Inostroza
    orcid: 0000-0002-4708-3275
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Vinicius Schaedler Damin
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 2
  - name: Erida Gjini
    orcid: 0000-0002-0493-3102
    affiliation: 3
affiliations:
 - name: Department of Biomedical Engineering, Instituto Superior Tecnico, University of Lisbon, Lisbon, Portugal
   index: 1
 - name: Department of Electrical and Computer Engineering, Instituto Superior Tecnico, University of Lisbon, Lisbon, Portugal
   index: 2
 - name: Center for Computational and Stochastic Mathematics, Instituto Superior Tecnico, University of Lisbon, Lisbon, Portugal
   index: 3
date: 30 June 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
[NAME]: MetaSpread
---

# Summary

*Keywords:* cancer, growth, metastatic spread, multi-scale dynamics, simulation

# Statement of need

Metastasis is a complex biological process that varies greatly between individuals and oncological phenotype. Currently, there is no easy-to-use framework that allows the use of real patients data for *in silico* metastasis spread modelling. Thus, a free, accesible and extendable software package is needed for providing a stable framework for metatasis simulation.

# Cancer growth and spread model

We provide a simulation framework in Python for the mathematical model in [@Franssen2019]. This paper proposed a hybrid modeling framework where cellular growth and metastatic spread are described and simulated in a spatially explicit manner, accounting for deterministic and stochastic dynamics. The model incorporates several key processes such as the interaction between epithelial and mesenchymal cells, the role of the extracellular matrix, diffusion, haptotaxis, circulation of cancer cells in the vasculature and seeding and growth in secondary sites. A 2D model of cancer is developed in Python. For this, we consider two different phenotypes: the Mesenchymallike cancer cells, which are capable of rapidly diffusing through the tissue and can break through the vasculature, but they reproduce at a minor rate compared to the Epithelial-like cells, which in contrast reproduce at a higher rate, but diffuse much more slowly and also can not break through the vasculature wall alone. Their movement is modeled considering 2 partial differential equations. The model includes two more equations: one for the generation of MMP-2, a chemical that favours the spread of cancer cells, and another for the degradation of the extracellular
matrix, which also favours the haptotactic movement of the cancer cells. For the simulation of the spatiotemporal growth dynamics, the system of PDE’s is discretized, and several 2-dimensional grids are established, representing the primary site and the metastatic sites, in which only the primary site is seeded with an initial number and distribution of cells. In order for the cells to migrate to another site, they must travel through the vasculature, which they do if they intravasate by one of the several randomly selected points in the grid that represent entrances to the vasculature system (see Figure 1). The extravasation to one of the metastatic sites only occurs if they survive. The model is programmed using MESA, [@python-mesa-2020] a Python Package for Agent-based modeling.

# Structure of the simulation platform

## API

For user interaction, the simulation platform is composed on three main modules:

- **Parameters.py:** it contains all the parameters detailed in section **Siulation parameters.** The user can edit this file directly before running a new simulation.

- **Simulation.py:** the file used for generating a new simulation. It will read the parameters found in the file ``parameters.py``.

## Visualization and analysis

- **graphsGenerator.py:** this will generate the graphs for the cancer cell distribution in all of the simulation's sites, and the ECM and MMP concentrations. If no parameters are given, it will create the graphs for all of the available time entries.

- **videoGenerator.py:** it will generate a video showing how the simulation changed through time.

- **analysis.py:** it will generate different graphs for analysis after the simulation.

## Background structure

There are three main classes that control the model's behavior, each contained in their homonym ``.py`` file.

- **CancerModel.py:** The whole model instances from this class. It governs the initial conditions, and the calculation of the MMP2 and ECM concentrations on each step.

- **CancerCell:** Each individual cell is an instance of this class. The movement rules for the cells are controlled here, and the intravasation of the cells to the vasculature system when they go into contact with a vessel point.

- **Vessel:** The vessels are modeled as MESA agents too. They do not posses any property, they only objective is to occupy one grid space, so that when a cancer cell moves to its position, it intravasates.

### Flow diagram

The algorithm that will be executed on each step is detailed in the following flow diagram:

# Simulation parameters

# Outlook

Possible extensions:

- Computational challenges (interface with a webpage, etc)

- EMT TME transformation

- Evolution or mutations of cells

- Immune system

- Healthy tissue cells

- Interaction with drugs

# Download, installation and run details

# Mathematical model

The dimensionless model, as described by Franssen et al. [@Franssen2019] corresponds to:




# # Acknowledgements

We acknowledge contributions fromMMurillo Teixeira

# References