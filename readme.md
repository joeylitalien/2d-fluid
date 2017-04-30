# Poisson Solvers for 2D Fluid Simulation
This repository serves for the final project of Paul Kry's Physically Based Animation course at McGill University, Winter 2017.

## Numerical Solvers
* Gauss-Seidel relaxation
* Successive over-relaxation
* Conjugate gradient
* Preconditioned conjugate gradient with incomplete Cholesky factorization

## Comparison
To compare algorithms, first run Gauss-Seidel with a very high number of iterations to establish ground truth for 60 seconds. Then, run the same scene with any solver for the same amount of time. The mean square error is output in the console.

## Controls
|Control           | Action                               |
|:-----------------|:--------------------------------------|
| Mouse left-click | Add source                           |
| Mouse right-lick | Delete source                        |
| Arrow up/down    | Increase/decrease source temperature |
| Space            | Play/pause                           |
| Enter            | Record                               |
| S                | Increase time                        |
| R                | Reset scene                          |
| N                | Delete all sources                   |
| C                | Add filament                         |
| 1-2-3-4          | Build appropriate test cases         |
