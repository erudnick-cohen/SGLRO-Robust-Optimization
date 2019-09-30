Non-Convex Robust Opimization Algorthm: Scenario Generation with Local Robust Optimization (SGLRO)

Robust optimization is a type of optimization that accounts for uncertainty by finding a solution which is feasible under all possible uncertain parameter values. 
This project implements a solver for robust optimization problems which contain non-convex objective functions or constraints.


This is an implementation of the Scenario Generation with Local Robust Optimization (SGLRO) algorithm presented in Rudnick-Cohen et al. 2019.
SGLRO is a sampling based approach, it randomly samples scenarios and uses them to generate worst case scenarios, which allows it to find a robust optimal solution. It also uses a local robust optimization step to ensure its final solution is correct.

The function SGLRO runs the SGLRO algorithm, see the file SGLRO.m for a list of inputs to this function and what they do.

The examples folder contains code for all the examples from Rudnick-Cohen et al. 2019, which demonstrate how to use SGLRO.m.

Sample code is provided for several other robust optimization approaches if you wish to compare results against them. Note that the code for these implementations is a bit messy relative to the SGLRO implementation.

If you use this project in an academic work, please cite this project using its associated journal publication:

E. Rudnick-Cohen, J. W. Herrmann, and S. Azarm, “Non-convex feasibility robust optimization via scenario generation and local refinement” ASME Journal of Mechanical Design, 2019.

