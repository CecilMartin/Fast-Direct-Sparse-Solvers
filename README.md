# Fast Direct Sparse Solvers

### Members: [Zhe Chen](zc1291@cims.nyu.edu); [Guanchun Li](guanchun.li@nyu.edu)



Fast Solver Course Final Project, 2020 Fall, Courant, NYU

Chap. 21 & 22 of Fast Direct Solvers for Elliptic PDEs by Per-Gunnar Martinsson.

---

### Intro

Guanchun and I(Zhe Chen) implemeeted this matlab library for the FDSS(Fast direct sparse solvers) in 2D. Currently, Example for Laplace equation and Helmhotz equation are provides. And Fast solver with linear complexity such like sweeping scheme and buffered sweeping scheme are supported.


### Documents

1. [Final brief slides](https://github.com/CecilMartin/Fast-Direct-Sparse-Solvers/blob/master/doc/Fast_Solver_Final_Project_Pre__FDSS.pdf)

2. [report](ttps://github.com/CecilMartin/Fast-Direct-Sparse-Solvers/blob/master/doc/Fast_Solver_Final_Project.pdf)

3. [Intro slides for Chapter 21](https://github.com/CecilMartin/Fast-Direct-Sparse-Solvers/blob/master/doc/FDSS.pdf)

4. [Intro slides for Chapter 22](https://github.com/CecilMartin/Fast-Direct-Sparse-Solvers/blob/master/doc/Linear Complexity.pdf)


### Dependency

FLAM(Fast linear algebra in MATLAB: Algorithms for hierarchical matrices):

http://klho.github.io/FLAM/#algorithms

### Code structure

1. src/: Source codes such as sweeping scheme and buffered sweeping scheme.

2. example/: main scripts for running examples such as sweeping scheme and buffered sweeping scheme.

3. tools/:  auxilary functions.

4. main.m: example script of one FDSS run.

5. main_plot.m: scripts to generate the figures.

### Usage:

1. Select the function get_A (stiffness matrix), get_f (load vector) to satisfy your problem. bv_fun to provide the boundary value function handle.

2. Modified main.m or scripts in example/ to get a main script. Change the parameters such like n

3. Change the FLAM_path to specify the directory where you install FLAM and it will automatically set up the path environment.



