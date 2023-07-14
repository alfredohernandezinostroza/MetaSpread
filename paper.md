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
date: 30 June 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

*Keywords:* cancer, growth, metastatic spread, multi-scale dynamics, simulation

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Cancer growth and spread model

We provide a simulation framework in Python for the mathematical model in [@Franssen2019], ,. This paper proposed a hybrid modeling framework where cellular growth and metastatic spread are described and simulated in a spatially explicit manner, accounting for deterministic and stochastic dynamics. The model incorporates several
key processes such as the interaction between epithelial and mesenchymal cells, the role of the extracellular matrix, diffusion, haptotaxis, circulation of cancer cells in the vasculature and seeding and growth in secondary sites. A 2D model of cancer is developed in Python. For this, we consider two different phenotypes: the Mesenchymallike cancer cells, which are capable of rapidly diffusing through the tissue and can break through the vasculature, but they reproduce at a minor rate compared to the Epithelial-like cells, which in contrast reproduce at a higher rate, but diffuse much more slowly and also can not break through the vasculature wall alone. Their movement is modeled considering 2 partial differential equations. The model includes two more equations: one for the generation of MMP-2, a chemical that favours the spread of cancer cells, and another for the degradation of the extracellular
matrix, which also favours the haptotactic movement of the cancer cells. For the simulation of the spatiotemporal growth dynamics, the system of PDE’s is discretized, and several 2-dimensional grids are established, representing the primary site and the metastatic sites, in which only the primary site is seeded with an initial number and distribution of cells. In order for the cells to migrate to another site, they must travel through the vasculature, which they do if they intravasate by one of the several randomly selected points in the grid that represent entrances to the vasculature system (see Figure 1). The extravasation to one of the metastatic sites only occurs if they survive. The model is programmed using MESA, [@python-mesa-2020] a Python Package for Agent-based modeling.

# Structure of the simulation platform

### Flow diagram

The algorithm that will be executed on each step is detailed in the following flow diagram:

# Simulation parameters

$$
\begin{array}{|c|c|c|c|}
\hline
 & \text{Variable name} & \text{ Description } & \text{Value}  \\
\hline
\Delta t & \texttt{th} & \text{ Time step } & 1\times 10^{-3}  \\
\hline
\begin{array}{ c }
\Delta x\\
\Delta y
\end{array} & \begin{array}{ c }
\texttt{xh}\\
\texttt{yh}
\end{array} & \text{ Space step } & 5\times 10^{-3}  \\
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
\Phi _{\mathrm{E}} & \texttt{phiE} & \text{ Epithelial hapto tactic sensi tivity coefficient } & 5\times 10^{-4}  \\
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
T_{V} & \texttt{vasculature\_time} & \text{ Time CTCs spend in the vasculature } & 0.18  \\
\hline
T_{\mathrm{M}} & \texttt{doublingTimeE}& \text{ Epithelial doubling time } & 3  \\
\hline
T_{\mathrm{E}} & \texttt{doublingTimeM} & \text{ Mesenchymal doubling time } & 2  \\
\hline
\mathcal{P}_{s} & \texttt{single\_cell\_survival} & \text{ Single CTC survival probability } & 5\times 10^{-4}  \\
\hline
\mathcal{P}_{C} & \texttt{cluster\_survival} & \text{ CTC cluster sunvival probability } & 2.5\times 10^{-2}  \\
\hline
\mathcal{E}_{1} & \texttt{E1} & \text{ Extravasation probability to bones } & \sim 0.5461  \\
\hline
\mathcal{E}_{2} & \texttt{E2} & \text{ Extravasation probability to lungs } & \sim 0.2553  \\
\hline
\mathcal{E}_{3} & \texttt{E3} & \text{ Extravasation probability to liver } & \sim 0.1986  \\
\hline
\end{array}
$$


# Simulation output, visualization and analysis

![Initial conditions of a sample run.\label{example-image-1}](example-image-1.png)

![Example results for the amount of cells in the vasculature\label{example-image-2}](example-image-2.png)

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

$$
\begin{aligned}
\frac{\partial c_{\mathrm{E}}}{\partial t} & =D_{\mathrm{E}} \nabla ^{2} c_{\mathrm{E}} -\Phi _{\mathrm{E}} \nabla \cdot ( c_{\mathrm{E}} \nabla w)\\
\frac{\partial c_{\mathrm{M}}}{\partial t} & =D_{\mathrm{M}} \nabla ^{2} c_{\mathrm{M}} -\Phi _{\mathrm{M}} \nabla \cdot ( c_{\mathrm{M}} \nabla w)\\
\frac{\partial m}{\partial t} & =D_{m} \nabla ^{2} m+\Theta c_{\mathrm{M}} -\Lambda m\\
\frac{\partial w}{\partial t} & =-( \Gamma _{1} c_{\mathrm{M}} +\Gamma _{2} m) w
\end{aligned}
$$

Where, once discretized take up the form:

$$
\begin{aligned}
c_{Ei,j}^{n+1} = & \mathcal{P}_{0} c^{n}_{Ei-1,j} +\mathcal{P}_{1} c^{n}_{Ei+1,j} +\mathcal{P}_{2} c^{n}_{Ei,j+1} +\mathcal{P}_{3} c^{n}_{Ei,j-1} +\mathcal{P}_{4} c^{n}_{Ei,j}\\
c_{Mi,j}^{n+1} = & \mathcal{P}_{0} c^{n}_{Mi-1,j} +\mathcal{P}_{1} c^{n}_{Mi+1,j} +\mathcal{P}_{2} c^{n}_{Mi,j+1} +\mathcal{P}_{3} c^{n}_{Mi,j-1} +\mathcal{P}_{4} c^{n}_{Mi,j}\\
m_{i,j}^{n+1} = & D_{m}\frac{\Delta t_{a}}{( \Delta x_{a})^{2}}\left( m_{i+1,j}^{n} +m_{i-1,j}^{n} +m_{i,j+1}^{n} +m_{i,j-1}^{n}\right)\\
 & +m_{i,j}^{n}\left( 1-4D_{m}\frac{\Delta t_{a}}{( \Delta x_{a})^{2}} -\Delta t\Lambda \right) +\Delta t_{a} \Theta c^{n}_{Mi,j}\\
w_{i,j}^{n+1} = & w_{i,j}^{n}\left[ 1-\Delta t_{a}\left( \Gamma _{1} c{_{M}^{n}}_{i,j} +\Gamma _{2} m_{i,j}^{n}\right)\right]
\end{aligned}
$$

Where:

$$
\begin{aligned}
\mathcal{P}_{0} : & \mathcal{P}_{i-1,j}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} -\frac{\Phi _{k}}{4}\left( w_{i+1,j}^{n} -w_{i-1,j}^{n}\right)\right]\\
\mathcal{P}_{1} : & \mathcal{P}_{i+1,j}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} +\frac{\Phi _{k}}{4}\left( w_{i+1,j}^{n} -w_{i-1,j}^{n}\right)\right]\\
\mathcal{P}_{2} : & \mathcal{P}_{i,j+1}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} +\frac{\Phi _{k}}{4}\left( w_{i,j+1}^{n} -w_{i,j-1}^{n}\right)\right]\\
\mathcal{P}_{3} : & \mathcal{P}_{i,j-1}^{n} :=\frac{\Delta t}{(\Delta x)^{2}}\left[ D_{k} -\frac{\Phi _{k}}{4}\left( w_{i,j+1}^{n} -w_{i,j-1}^{n}\right)\right]\\
\mathcal{P}_{4} : & \mathcal{P}_{i,j}^{n} :=1-(\mathcal{P}_{0} +\mathcal{P}_{1} +\mathcal{P}_{2} +\mathcal{P}_{3})
\end{aligned}
$$
# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:

- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References