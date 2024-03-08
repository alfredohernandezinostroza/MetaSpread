---
title: 'MetaSpread: A cancer growth and metastatic spread simulation program in Python'
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
    equal-contrib: true
    affiliation: 1
  - name: Vinicius Schaedler Damin
    orcid: 0000-0000-0000-0000
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Erida Gjini
    affiliation: 3
    corresponding: true # (This is how to denote the corresponding author)
affiliations:
 - name: Department of Biomedical Engineering, Instituto Superior Tecnico, University of Lisbon, Lisbon, Portugal
   index: 1
 - name: Department of Electrical and Computer Engineering, Instituto Superior Tecnico, University of Lisbon, Lisbon, Portugal
   index: 2
 - name: Center for Computational and Stochastic Mathematics, Instituto Superior Tecnico, University of Lisbon, Lisbon, Portugal
   index: 3
date: 9 March 2023
bibliography: paper.bib

---

# Summary

We develop and provide MetaSpread, an open source simulation framework and interactive program in Python, based on a mathematical model of tumor growth and metastatic spread, designed by [franssen2019]. This paper proposed a hybrid modeling and computational framework where cellular growth and metastatic spread are described and simulated in a spatially explicit manner, accounting for stochastic individual cell dynamics and deterministic dynamics of abiotic factors. This model incorporates several key processes such as the growth and movement of epithelial and mesenchymal cells, the role of the extracellular matrix, diffusion, haptotaxis, circulation and survival of cancer cells in the vasculature, and seeding and growth in secondary sites. In the software that we develop, these growth and metastatic dynamics are programmed using MESA, @[python-mesa-2020] a Python Package for Agent-based modeling.

*Keywords:* cancer, growth, metastatic spread, multi-scale dynamics, simulation

# Statement of need

Models of tumor growth and metastatic spread are key for understanding the key underlying biological processes and clinical evolution in patients. Mathematical models can be of different level of detail, computational or theoretical, spatial or non-spatial in nature, and can have several mechanisms explicit or implicitly embedded in them, including interaction with resources, biomechanical signals, cellular competition, mutation and migration @[waclaw_spatial_2015; franssen2021novel;macnamara2020computational;chaplain2020multiscale;sadhukhan2022multi;opasic2020cancersim]. While theoretical and analytical advances remain crucial in mathematical models of cancer, computational approaches that offer direct simulation platforms for efficient numerical exploration, focused study and hypothesis testing are also very much needed. Here, we contribute to this aspect, by offering an open source simulation framework in Python for spatio-temporal progression of tumor and metastatic spread. We build the simulation framework on a hybrid mathematical model developed by @[franssen2019], extending the previous work by @[anderson1998continuous; anderson2000mathematical], in close agreement with empirical data @[sabeh2009protease;newton2015spatiotemporal]. This contribution aims to bridge the gap between mathematicians, oncologists, biologists, computer scientists and interested researchers working in the field of cancer metastatic progression, in need of a computational framework for interdisciplinary study and collaboration.

# Cancer growth and spread model

A 2-dimensional multigrid hybrid spatial model of cancer dynamics is developed in Python. Here we combine the stochastic individual based dynamics of singles cells with deterministic dynamics of the abiotic factors. The algorithm for dynamic progression at each time step is depicted in Figure . In the tumor site we consider two different cancer cell phenotypes: epithelial (epithelial-like) and mesenchymal (mesenchymal-like) cells. The epithelial-like cancer cells reproduce at a higher rate, but diffuse more slowly than mesenchymal cells, which reproduce at a lower rate but diffuse more rapidly. Furthermore, epithelial cells cannot break through the vasculature wall alone, as they require the presence of mesenchymal cells to be able to intravasate. The cellular growth and movement in space is modeled considering 2 partial differential equations, where random (diffusion) and non-random (haptotaxis) movement are implemented. The model includes two additional equations: one for the spatio-temporal dynamics of matrix metalloproteinase 2 (MMP-2), a chemical that favors the spread of cancer cells, and another for the degradation of the extracellular matrix (ECM), which also favors the haptotactic movement of the cancer cells. 
For the simulation of the spatio-temporal growth dynamics, and metastatic spread, the system of PDE's is discretized, and several 2-dimensional grids are established, representing the primary site and the metastatic sites. Only the primary site is seeded with an initial number and distribution of cells. In order for the cells to migrate to another site, they must travel through the vasculature, which they do if they intravasate by one of the several randomly selected points in the grid that represent entrances to the vasculature system. The extravasation to one of the metastatic sites only occurs if they survive, a process that is modeled with net probabilistic rules considering time spent in the vasculature, cluster disaggregation, cell type, and potential biases to different destinations. 
The dimensionless model, as described by @[franssen2019] in Appendix A of their paper, corresponds to 4 PDEs, where the key variables reflect local densities of epithelial cells ($c_E$) and mesenchymal cells ($c_M$), and concentrations of MMP2 ($m$) and extracellular matrix ($w$):

$$
\begin{aligned}
\frac{\partial c_{\mathrm{E}}}{\partial t} & =D_{\mathrm{E}} \nabla ^{2} c_{\mathrm{E}} -\Phi _{\mathrm{E}} \nabla \cdot ( c_{\mathrm{E}} \nabla w)\\
\frac{\partial c_{\mathrm{M}}}{\partial t} & =D_{\mathrm{M}} \nabla ^{2} c_{\mathrm{M}} -\Phi _{\mathrm{M}} \nabla \cdot ( c_{\mathrm{M}} \nabla w)\\
\frac{\partial m}{\partial t} & =D_{m} \nabla ^{2} m+\Theta c_{\mathrm{M}} -\Lambda m\\
\frac{\partial w}{\partial t} & =-( \Gamma _{1} c_{\mathrm{M}} +\Gamma _{2} m) w
\end{aligned}
$$

Discretizing equations 1 and 2 in space and time, we obtain:

$$
\begin{aligned}
c_{Ei,j}^{n+1} = & \mathcal{P}_{0} c^{n}_{Ei-1,j} +\mathcal{P}_{1} c^{n}_{Ei+1,j} +\mathcal{P}_{2} c^{n}_{Ei,j+1} +\mathcal{P}_{3} c^{n}_{Ei,j-1} +\mathcal{P}_{4} c^{n}_{Ei,j}\\
c_{Mi,j}^{n+1} = & \mathcal{P}_{0} c^{n}_{Mi-1,j} +\mathcal{P}_{1} c^{n}_{Mi+1,j} +\mathcal{P}_{2} c^{n}_{Mi,j+1} +\mathcal{P}_{3} c^{n}_{Mi,j-1} +\mathcal{P}_{4} c^{n}_{Mi,j}\\
\end{aligned}
$$

Where $\mathcal{P}_0$ to $\mathcal{P}_4$:

$$
\begin{aligned}
\mathcal{P}_{0} : & \mathcal{P}_{i-1,j}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} -\frac{\Phi _{k}}{4}\left( w_{i+1,j}^{n} -w_{i-1,j}^{n}\right)\right]\\
\mathcal{P}_{1} : & \mathcal{P}_{i+1,j}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} +\frac{\Phi _{k}}{4}\left( w_{i+1,j}^{n} -w_{i-1,j}^{n}\right)\right]\\
\mathcal{P}_{2} : & \mathcal{P}_{i,j+1}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} +\frac{\Phi _{k}}{4}\left( w_{i,j+1}^{n} -w_{i,j-1}^{n}\right)\right]\\
\mathcal{P}_{3} : & \mathcal{P}_{i,j-1}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} -\frac{\Phi _{k}}{4}\left( w_{i,j+1}^{n} -w_{i,j-1}^{n}\right)\right]\\
\mathcal{P}_{4} : & \mathcal{P}_{i,j}^{n} :=1-(\mathcal{P}_{0} +\mathcal{P}_{1} +\mathcal{P}_{2} +\mathcal{P}_{3})
\end{aligned}
$$

represent the probabilities for a cell to move up, down, left, right, or stay in place, and where $k=E,M$ can refer to an epithelial-like or mesenchymal-like cell. Each cell on every grid point $(x_i,y_j)$ is modeled as an individual agent, which obeys probability rules for growth and movement. There is a maximal carrying capacity for each grid point given by $Q,$ (assumed equal to 4 in @[franssen2019]), to represent competition for space. There exist a doubling time $T_E$ and $T_M$ for epithelial and mesenchymal cells at which all the cells present in all grids will reproduce, duplicating in place, but never exceeding $Q$.

For the abiotic factors the discretization takes the form (see Appendices in @[franssen2019]):

$$
\begin{aligned}
m_{i,j}^{n+1} = & D_{m}\frac{\Delta t_{a}}{( \Delta x_{a})^{2}}\left( m_{i+1,j}^{n} +m_{i-1,j}^{n} +m_{i,j+1}^{n} +m_{i,j-1}^{n}\right)\\
 & +m_{i,j}^{n}\left( 1-4D_{m}\frac{\Delta t_{a}}{( \Delta x_{a})^{2}} -\Delta t\Lambda \right) +\Delta t_{a} \Theta c^{n}_{Mi,j}\\
w_{i,j}^{n+1} = & w_{i,j}^{n}\left[ 1-\Delta t_{a}\left( \Gamma _{1} c{_{M}^{n}}_{i,j} +\Gamma _{2} m_{i,j}^{n}\right)\right]
\end{aligned}
$$

Where $i,j$ reflect the grid ($i,j$) and $n$ the time-point. In this discretization two different time and spatial steps are used for the cell population and the abiotic factors (ECM and MMP-2), namely delta t and delta x, and delta t a and delta x a, respectively.

# Simulation parameters

The biological parameters of the model and the simulation values are summarized in Table 1, tailored to breast cancer progression and early-stage dynamics prior to any treatment and in a pre-angiogenic phase (less than 0.2 cm in diameter). We provide the default values used by @[franssen2019], as informed by biological and empirical considerations (see also Table 1 and references therein in @[franssen2019franssen2019mics represent a two-dimensional cross-section of a small avascular tumor and run on a 2-dimensional discrete grid (spatial domain $[0,1] \times [0,1]$ corresponding to physical domain of size $[0,0.2]\text{ cm} \times [0,0.2]\text{ cm}$), where each grid element corresponds to a spatial unit of dimension $(\Delta x,\Delta y)$, and where position $x_i,y_j$ corresponds to $i \Delta x$ and $j \Delta y$. Cancer cells are modeled as discrete agents whose growth and migration dynamics follow probabilistic rules, whereas the abiotic factors MMP2 and extracellular matrix dynamics follow the deterministic PDE evolution, discretized by an explicit five-point central difference discretization scheme together with zero-flux boundary conditions. The challenge of the simulation lies in coupling deterministic and agent-based stochastic dynamics, and in formulating the interface between the primary tumor Grid 1 and the metastatic sites (Grids 2,..$k$). Each grid shares the same parameters, but there can be biases in connectivity parameters between grids ($\mathcal{E}_{k}$ parameters).

Cell proliferation is implemented in place by generating a new cell when the doubling time is completed, for each cell in each grid point. But if the carrying capacity gets surpassed, then there is no generation of a new cell. The movement of the cells is implemented through the probabilities in Equations, which are computed at each time point and for each cell and contain the contribution of the random diffusion process and non-random haptotactic movement. If a cell lands in a grid point that contains a vasculature entry point, it is typically removed from the main grid and added to the vasculature. But there are details regarding the type of cells (E or M) and vasculature entry points (normal or ruptured) further described by @[franssen2019].

The vasculature is the structure connecting the primary and secondary sites, and it represents a separate compartment in the simulation framework. Single cells or clusters of cells can enter the vasculature either through a ruptured or normal vessel, and they can remain there for a fixed number of time $T_V$, representing the average time a cancer cell spends in the blood system. Each cell belonging to a cluster in the vasculature can disaggregate with some probability. At the end of the residence time in the vasculature, each cell's survival is determined randomly with probabilities that are different for single and cluster cells, and the surviving cells are randomly distributed on the secondary sites. To implement this vasculature dynamics in the algorithm, the vasculature is represented as a dictionary where the keys refer to the time-step in which there are clusters ready to extravasate. Intravasation at time $t$ corresponds to saving the cells into the dictionary with the associated exit time $t+T_V$. 

Extravasation rules follow the setup in the original paper @[franssen2019], ensuring arriving cells do not violate the carrying capacity. Metastatic growth after extravasation follows the same rules as in the original grid. 