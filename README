Global Model — Multicomponent Pool Fires
=========================================

MATLAB implementation of the transient "global model" for multicomponent
liquid pool fires described in:

Yip, A., Haelssig, J.B., Pegg, M.J. (2021). "Multicomponent pool fires:
Trends in burning rate, flame height, and flame temperature." Fuel, 284,
118913. https://doi.org/10.1016/j.fuel.2020.118913

Background
----------
Multicomponent fuel mixtures (e.g. ethanol-water, ethanol-isopropanol,
ethanol-hexane) undergo preferential distillation of their more volatile
components as a pool fire burns, so the liquid composition, and therefore
the burning rate, flame height, flame temperature, and soot production,
all evolve over time. The model tracks these coupled changes by treating
the pool as a well-mixed liquid in vapour-liquid equilibrium (via the
Wilson activity-coefficient model) and solving a system of ODEs for the
mass, composition, and flame temperature of the fire as it depletes.
Empirical correlations (Babrauskas, Ditch et al., Heskestad) are used
alongside the liquid-phase model to predict burning rate, flame height,
and flame temperature, and results are compared against the fire
experiments reported in the paper.

Repository layout
------------------
- `Global Model Multicomponent/` — main entry point. `main.m` sets up the
  fuel mixture, initial conditions, and solves the transient global model;
  `run.m` runs `main.m` and plots/tabulates the results (burning rate,
  flame temperature, mass, composition over time). Currently supports
  ethanol-water and ethanol-heptane mixtures. Also includes VLE test
  scripts (`testVLE_*.m`) and mass/volume fraction conversion utilities.
- `Library/` — supporting classes and functions used by the model:
  - `Species/` — `Species`/`Mixture` classes and property definitions for
    the pure components and air used in the model (ethanol, water,
    heptane, hexane, propanol, isopropanol, oxygen, nitrogen, air,
    carbon dioxide).
  - `Functions/` — calculations for flame height, flame temperature,
    radiation, mass/mole flow rates, and VLE data.
- `Tests/` — unit test scripts (`test_*.m`) that verify individual
  functions and utilities (e.g. flame height, flame temperature, burning
  rate, convection, air entrainment) used by the global model.

Usage
-----
1. In `Global Model Multicomponent/main.m`, specify the fuel mixture,
   initial mass fractions, and initial mass. Only the following binary
   mixtures are currently supported (species must be listed in this
   order): (i) Ethanol-Water, (ii) Ethanol-Heptane.
2. Run `Global Model Multicomponent/run.m` to solve the model and generate
   plots of burning rate, flame temperature, and composition over time.
