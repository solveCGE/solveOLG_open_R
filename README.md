# solveOLG open economy in R
Solves a simple AK-OLG-model for a small open economy in R

## About
Shows how to solve a simple deterministic overlapping-generations model (OLG) of Auerbauch-Kotlikoff type, solving for the transition path between two steady-states. The code is not optimized for speed, makes excessive use of global variables, and is mainly meant for instructional purposes. However, the algorithm as such is efficient and written with parallelization in mind (at every iteration the household problems of all finitely lived representative cohorts could be computed in parallel). For larger models, it is recommended to use a more efficient language that supports shared memory multithreading and pass-by-reference, e.g. C++ with the Armadillo library or Julia.

A model description can be found here: <https://github.com/solveCGE/solveOLG_doc>.

## How to run
Parameters can be set in `calib.R`. Policy shocks are defined in `run.R` (or just uncomment some of the predefined exemplary shocks). The model is then solved by just running `run.R`. 

## Author
Philip Schuster
