# PoissonSolver

## Introduction
<img src="results/vecteurPropreLogo.png" width="200" align="right">
This project presents a numerical solution to the two-dimensional Poisson equation using a C program. The solution involves defining arbitrary membrane shapes, discretizing the domain, and employing a Multi-Grid method for iterative problem-solving. Further, the project explores the optimization of the algorithm through a relaxation parameter and enhances solver performance using a Multi-Grid preconditioner with the PRIMME solver.

## Program Structure
The program is structured into 20 files, each dedicated to specific functions. The main file, *main.c*, offers four distinct options:

- **The first option** includes a *two-grid method* originally designed for a particular problem, though the second and third options can also implement this method.
- **The second option** enables the use of a *multi-grid solver* that can operate in either V or W cycles, with a variable number of levels.
- **The third option** establishes a *multi-grid preconditioner* that, at the end of each iteration, applies a correction to the solution, incorporating a *relaxation factor*, and calculates both the minimum/maximum eigenvalues and the optimal relaxation factor. It also checks the algorithm's direct stability.
- **The fourth option** employs the *PRIMME solver* to tackle an eigenvalue problem, with the *multi-grid preconditioner* enhancing the solver's convergence rate.

## Discretization
The project involves solving the Poisson equation for a given membrane $\Omega \subset \mathbb{R}^2$, with Dirichlet boundary conditions specified on $\partial \Omega$.
$$
    \begin{cases} 
        \begin{aligned}
            -\Delta u &= 0 &&\text{on} \Omega \subset \mathbb{R}^2\\ 
            u &= e^{\sqrt{x^2+y^2}} &&\text{on} \partial \Omega\\
        \end{aligned}
    \end{cases}
$$


## Multi-Grid Method
<figure>
  <img src="results/P8_twoGrid_differentSmoothing.png" width="400" align="center">
  <figcaption>Convergence of the Two-Grid Method for Residuals Based on the Number of Pre- and Post-Smoothing Iterations (641 Discretization Points)
</figure>
The Multi-Grid method enhances solution convergence through pre-/post-smoothing iterations and corrections on coarser grids. The method's effectiveness and the impact of various smoothing iterations on convergence are demonstrated.

## Preconditioning with Multi-Grid
The use of a Multi-Grid preconditioner with the PRIMME solver for eigenvalue problems showcases significant improvements in solver performance.


## Conclusion
The project effectively demonstrates the application of the Multi-Grid method to solve the 2D Poisson equation and optimize eigenvalue problem solving with the PRIMME solver.
