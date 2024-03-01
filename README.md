# PoissonSolver
<p float="left">
  <img src="/results/solFinaleCercle.png" width="300" align="center" />
  <img src="/results/residuMultiGridP8_m641.png" width="300" align="center" /> 
  <img src="/results/solutionP15.png" width="300" align="center"/>
</p>

## Introduction
This project presents a numerical solution to the two-dimensional Poisson equation using a C program. The solution involves defining arbitrary membrane shapes, discretizing the domain, and employing a Multi-Grid method for iterative problem-solving. Further, the project explores the optimization of the algorithm through a relaxation parameter and enhances solver performance using a Multi-Grid preconditioner with the PRIMME solver.

## Program Structure
The program comprises 20 files, each serving specific functions. The main file, \textit{main.c}, offers four options ranging from a two-grid method to employing a Multi-Grid preconditioner for the PRIMME solver.

## Discretization
The project involves solving the Poisson equation for a given membrane $\Omega \subset \mathbb{R}^2$, with boundary conditions specified on $\partial \Omega$. The domain is discretized, and the solution $u$ at each grid point is approximated using the finite difference method, leading to a linear system $Au = b$.

## Multi-Grid Method
The Multi-Grid method enhances solution convergence through pre-/post-smoothing iterations and corrections on coarser grids. The method's effectiveness and the impact of various smoothing iterations on convergence are demonstrated.

## Preconditioning with Multi-Grid
The use of a Multi-Grid preconditioner with the PRIMME solver for eigenvalue problems showcases significant improvements in solver performance.

<img src="results/discretisation_domaine.png" width="400" align="center">

<img src="results/solutionP15" width="400" align="center">


## Conclusion
The project effectively demonstrates the application of the Multi-Grid method to solve the 2D Poisson equation and optimize eigenvalue problem solving with the PRIMME solver.
