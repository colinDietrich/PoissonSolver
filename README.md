# PoissonSolver

## Introduction
<img src="results/vecteurPropreLogo.png" width="200" align="right">
This project presents a numerical solution to the two-dimensional Poisson equation using a C program. The solution involves defining arbitrary membrane shapes, discretizing the domain, and employing a Multi-Grid method for iterative problem-solving. Further, the project explores the optimization of the algorithm through a relaxation parameter and enhances solver performance using a Multi-Grid preconditioner with the <a href="https://www.cs.wm.edu/~andreas/software/">PRIMME</a> solver.

## Program Structure
The program is structured into 20 files, each dedicated to specific functions. The main file, `main.c`, offers four distinct options:

- **The first option** includes a *two-grid method* originally designed for a particular problem, though the second and third options can also implement this method.
- **The second option** enables the use of a *multi-grid solver* that can operate in either V or W cycles, with a variable number of levels.
- **The third option** establishes a *multi-grid preconditioner* that, at the end of each iteration, applies a correction to the solution, incorporating a *relaxation factor*, and calculates both the minimum/maximum eigenvalues and the optimal relaxation factor. It also checks the algorithm's direct stability.
- **The fourth option** employs the *PRIMME solver* to tackle an eigenvalue problem, with the *multi-grid preconditioner* enhancing the solver's convergence rate.

## Discretization
The project involves solving the Poisson equation for a given membrane $\Omega \subset \mathbb{R}^2$, with Dirichlet boundary conditions specified on $\partial \Omega$.

$$
\begin{cases}
  \begin{aligned} 
    -\Delta u &= 0 &&\text{on} &&&\Omega \subset \mathbb{R}^2\\ 
    u &= e^{\sqrt{x^2+y^2}} &&\text{on} &&&\partial \Omega\\ 
  \end{aligned}
\end{cases}
$$

The domain is discretized with a discretization step $h = \frac{L}{m-1}$ where $L$ is the length of one side of the square membrane and $m$ is the number of points aligned in one direction of the grid. Thus, each of the grid points is defined as follows:

$$
\begin{aligned}
    (x_i,y_j) &= (ih,jh) &&\text{with} &&&i,j = 1,...,m-2
\end{aligned}
$$

The approximation of the solution $u$ is calculated at each of the points belonging to the interior of the domain $\Omega$ according to the finite difference:

$$
    \frac{4u_{i,j} - u_{i+1,j} - u_{i-1,j} - u_{i,j+1} - u_{i,j-1}}{h^2} = 0 + \mathcal{O}(h^2)
$$

where the edges of the domain are not taken into account because $u$ is already known there due to Dirichlet boundary conditions.
By gathering all the equations, we obtain the following linear system:

$$
    Au = b
$$

where $b$ is a null vector to which Dirichlet boundary conditions have been added.

## Multi-Grid Method Algorithm

```pseudo
FOR i = 0 to N
    - Apply ν1 pre-smoothing passes on A1u = b -> u1
    - Reduce the residue r2 = R1(b - A1u1)
        - Apply ν1 pre-smoothing passes on A2c = r2 -> c2
        - Reduce the residue r3 = R2(r2 - A2c2)
            - Apply ν1 pre-smoothing passes on A3c = r3 -> c3
            - Reduce the residue r4 = R3(r3 - A3c3)
                - Solve on the coarse grid: c4 = A4^(-1)r4
            - Prolong the correction: c3_tilde = c3 + P3c4
            - Apply ν2 post-smoothing passes on A3c = r3 -> c3'
        - Prolong the correction: c2_tilde = c2 + P2c3'
        - Apply ν2 post-smoothing passes on A2c = r2 -> c2'
    - Prolong the correction: u1_tilde = u1 + P1c2'
    - Apply ν2 post-smoothing passes on A1u = b -> u1'
END FOR

