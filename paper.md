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

Models of tumor growth and metastatic spread are key for understanding the key underlying biological processes and clinical evolution in patients. Mathematical models can be of different level of detail, computational or theoretical, spatial or non-spatial in nature, and can have several mechanisms explicit or implicitly embedded in them, including interaction with resources, biomechanical signals, cellular competition, mutation and migration @[waclaw_spatial_2015, franssen2021novel,macnamara2020computational,chaplain2020multiscale,sadhukhan2022multi,opasic2020cancersim]. While theoretical and analytical advances remain crucial in mathematical models of cancer, computational approaches that offer direct simulation platforms for efficient numerical exploration, focused study and hypothesis testing are also very much needed. Here, we contribute to this aspect, by offering an open source simulation framework in Python for spatio-temporal progression of tumor and metastatic spread. We build the simulation framework on a hybrid mathematical model developed by @[franssen2019], extending the previous work by @[anderson1998continuous, anderson2000mathematical], in close agreement with empirical data @[sabeh2009protease,newton2015spatiotemporal]. This contribution aims to bridge the gap between mathematicians, oncologists, biologists, computer scientists and interested researchers working in the field of cancer metastatic progression, in need of a computational framework for interdisciplinary study and collaboration.

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

# Structure of the simulation platform

The program can be run both interactively through the command line, or with explicit user command line arguments.

When run interactively, starting from the main menu, the following possibilities are offered: 

- **Run a new simulation:** the user can choose the *New Simulation* option to run a new simulation, with the arguments to be specified by the user being the maximal time for the dynamics, and the frequency of saving data (temporal resolution). Any other simulation parameter (see Table ) will be taken from the *simulation\_configs.csv* file in the main folder. At the end of the simulation the dynamics of the grids, including agents (cells and vasculature points), the vasculature dynamics and the MMP2 and ECM are saved in a properly identified directory, including a *configs.csv* recording the used parameters for this particular simulation. The file *CellsData.csv* in this directory will include all the information of all cells and vasculature points in the simulation, for every time step.

- **Load an existing simulation** The user can select *Load Simulation* from the main menu, and an existing simulation will be loaded, and continued for further time steps with the same parameters in its *configs.csv* file. The only parameters that the user has to select are the new temporal resolution and the maximum extra steps for the simulation to run.

- **Post-process data from a simulation** The generated *CellsData.csv* contains the information of every cancer cell at every time step and every grid of the simulation. In order to facilitate the study of the results, we provide with several post-processing methods: Data analysis, Graphical analysis and Video generation.

- **Data analysis:** several results will be summarized in *.csv* files, such as the vasculature and tumor dynamics.

- The files that account for cell growth, vasculature evolution, and tumor radius and diameter evolution consist of columns that register the state of a metric in each time step along the simulation. These easily allows plotting graphs of dynamics later on.

- The tumor growth files consist of 8 rows, where each 4 pairs correspond to the x and y coordinates of an agent in the simulation. The order of the rows from top to bottom correspond to mesenchymal, epithelial, normal vasculature points and ruptured vasculature points. These allow for easily plotting the positions of the agents, and thus, the state of the tumor, at each time step.

- The histogram files consists of two columns: one for the bins, and one for the frequency. The bins are for accounting the amount of grid points that have X amount of cells.

- **Graphical analysis:** in order to run this step, it is necessary to run the data analysis option first. When selected, plots of the tumor distribution, ECM, MMP-2 for each grid for all the amount of time points specified by the user will be produced. Furthermore, it will also produce other plots such as the dynamics of the cells in the vasculature, histograms of the cell number distribution over grid points, radius and diameter of the tumor over time, and cells number growth over time (total size of the tumor in each grid).

- **Video generation:** The user can choose the Videos option to gene~~~~rate videos animated from the consecutive snapshots of the spatial dynamics generated from the graphical analysis step.

# Simulation parameters

$$
\begin{array}{|c|c|c|c|}
\hline
 & \text{Variable name} & \text{ Description } & \text{Value}  \\
\hline
\Delta t & \texttt{th} & \text{ Time step } & 1\times 10^{-3}  \\
\hline
\Delta t_a & \texttt{tha} & \text{ Abiotic time step } & 1\times 10^{-3}  \\
\hline
\begin{array}{ c }
\Delta x\\
\Delta x_a
\end{array} & \begin{array}{ c }
\texttt{xh}\\
\texttt{xah}
\end{array} & \text{ Space step \\ Abiotic space step } & 5\times 10^{-3}  \\
\hline
D_{\mathrm{M}} & \texttt{dM} & \begin{array}{ c }
\text{ Mesenchymal-like cancer cell diffusion }\\
\text{ coefficient }
\end{array} & 1\times 10^{-4}  \\
\hline
D_{\mathrm{E}}& \texttt{dE} & \text{ Epithelial-like cancer cell diffusion coefficient } & 5\times 10^{-5}  \\
\hline
\Phi _{M} & \texttt{phiM} & \text{ Mesenchymal hap to tactic sensitivity coefficient } & 5\times 10^{-4}  \\
\hline
\Phi _{\mathrm{E}} & \texttt{phiE} & \text{ Epithelial hapto tactic sensitivity coefficient } & 5\times 10^{-4}  \\
\hline
D_{m} & \texttt{dmmp} & \text{ MMP-2 diffusion coefficient } & 1\times 10^{-3}  \\
\hline
\Theta & \texttt{theta} & \text{ MMP-2 production rate } & 0.195  \\
\hline
\Lambda & \texttt{Lambda} & \text{ MMP-2 decay rate } & 0.1  \\
\hline
\Gamma _{1} & \texttt{gamma1} & \text{ ECM degradation rate by MT1-MMP } & 1  \\
\hline
\Gamma _{2}& \texttt{gamma2} & \text{ ECM degradation rate by MMP-2 } & 1  \\
\hline
T_{V} & \texttt{vasculature\_time} & \text{ Steps CTCs spend in the vasculature } & 180  \\
\hline
T_{\mathrm{M}} & \texttt{doublingTimeE}& \text{ Epithelial doubling time } & 3  \\
\hline
T_{\mathrm{E}} & \texttt{doublingTimeM} & \text{ Mesenchymal doubling time } & 2  \\
\hline
\mathcal{P}_{s} & \texttt{single\_cell\_survival} & \text{ Single CTC survival probability } & 5\times 10^{-4}  \\
\hline
\mathcal{P}_{C} & \texttt{cluster\_survival} & \text{ CTC cluster survival probability } & 2.5\times 10^{-2}  \\
\hline
\mathcal{E}_{1,...,n} & \texttt{E1} & \text{ Extravasation probabilities} & \sim [0.54615461, 0.2553, 0.1986]  \\
\hline
\mathcal{P}_{d} & \texttt{disaggregation\_prob} & \text{ Individual cancer cell dissagregation probability} & 0.5  \\
\hline
Q & \texttt{carrying\_capacity} & \text{ Maximum amount of cells per grid point} & 4  \\
\hline
U_P & \texttt{normal\_vessels\_primary} & \text{ Amount of normal vessels present on the primary grid} & 2  \\
\hline
V_P & \texttt{ruptured\_vessels\_primary} & \text{ Amount of ruptured vessels present on the primary grid} & 8  \\
\hline
U_{2,...,n} & \texttt{secondary\_sites\_vessels} & \text{ Amount of vessels present on the secondary sites} & [10, 10]  \\
\hline
- & \texttt{n\_center\_points\_for\_tumor} & \text{ Amount of center-most grid points where the primary cells are going to be seeded} & 97  \\
\hline
- & \texttt{n\_center\_points\_for\_vessels} & \text{ Amount of center-most grid points where the vessels will not be able to spawn} & 200  \\
\hline
- & \texttt{gridsize} & \text{Length in gridpoints of the grid's side} & 201  \\
\hline
- & \texttt{grids\_number} & \text{ Amount of grids, including the primary site} & 3  \\
\hline
- & \texttt{mesenchymal\_proportion} & \text{ Proportion of mesenchymal-like cells to be seeded on the primary site} & 0.6  \\
\hline
- & \texttt{epithelial\_proportion} & \text{ Proportion of epithelial-like cells to be seeded on the primary site} & 0.4  \\
\hline
\end{array}
$$

# Simulation output, visualization and analysis

To illustrate the performance and capability of MetaSpread, we provide some figures and visualization of the simulations output. In Figure 2 we show a snapshot of cell distribution after approximately 5 days of growth. In Figure 3 we show a later snapshot of our simulations for cancer cell spread and ECM and MMP2 evolution. In Figure 4 we show temporal dynamics of summary variables, e.g. total cell counts over time up to 12.5 days, possible to be computed after simulation data post-processing. In movies S1-S2 we show how the simulation platform can be used for studying the biological effect of different perturbations in parameters. These movies illustrate animations of the spatiotemporal evolution of a tumor on the primary site in two cases: (S1) diffusion-dominated and (S2) haptotaxis-dominated cellular movement. The first leads to a regular spatiotemporal pattern of growth, more isotropic and round, the second leads to a more irregular growth over space with cellular protrusions extending in some directions.

![Initial conditions of a sample run.\label{example-image-1}](example-image-1.png)

<!-- ![Example results for the amount of cells in the vasculature\label{example-image-2}](example-image-2.png) -->

# Outlook

While the model originating from our program @[franssen2019] is simpler than later models developed for cancer invasion @[franssen2021novel,macnamara2020computational,chaplain2020multiscale], we believe the simple framework enables already deep study of the basic population dynamic processes involved in early tumor dynamics and metastatic growth, and engagement with interesting and important biology (reviewed in @[franssen2019franssen2019ient but not too hard level of complexity makes it a perfect tool for interaction by non-specialists in the mathematical field, medical doctors and for researchers willing to explore hypotheses with it, perform simulations or extract from it pedagogical value for students and the wider public. There are several directions for extensions of the algorithm and simulation package. These include improving the computational efficiency and speed of the simulation, which now requires about 24 hours for 28.000 time steps, corresponding to about 12 days. Another venue for extension could be including interaction with the immune system, developing explicitly the interaction of the cancer cells with healthy cells, implementing the effect of treatment, for example adaptive therapies @[west2023survey], mutations and EMT transition. Regarding the metastatic spread, a novelty would be to consider different parameters in different grids, allowing for differential suitability for growth and colonization by arriving cancer cells, which in the current formulation is captured only by the biases in arrival probabilities @[newton2015spatiotemporal]. On the computational side the main challenge relies on making the code flexible for parallel computing so that both the spatial and temporal resolution can be increased, and the scope of the phenomena investigated can be expanded, including cell-level heterogeneity. Finally, an interesting improvement would be the addition of an interactive interface that allows for better visualization of the tumor at different time steps, allowing for zooming in and out, changing the color of the cells, and more.

- Computational challenges (interface with a webpage, etc)

- Microbes playing an important role in cancer prognosis and development.

- EMT TME transformation

- Evolution or mutations of cells

- Immune system

- Healthy tissue cells

- Interaction with drugs

# Supporting information

Movie S1: Example 1 of spatiotemporal evolution of tumor growth in the primary site (default parameters, diffusion-dominated movement). Movie S2: Example 2 of spatiotemporal evolution of tumor growth in the primary site (parameters with haptotaxis-dominated movement of cells).



# Acknowledgements

We acknowledge the contribution of 

# References