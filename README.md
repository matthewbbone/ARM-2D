### Agent-based Modeling for Modular Theory Comparison

This repository contains the code for a presentation given at the [2023 International Conference of Computational Social Science](https://laura.alessandretti.com/public/pdf_accepted/paper86.pdf). In this project, the agent-based opinion dynamic model presented in [Axelrod et al.](https://www.pnas.org/doi/abs/10.1073/pnas.2102139118) was reimplemented and modified to operate with varying theoretical assumptions. Namely, the attraction/repulsion and homophily mechanisms were modified to operate on an ideological or affective basi, resulting in four versions of the model. These modulated versions of the model were then fit to poll data from the American National Election Study (ANES) to see which was able to provide the closest approximation.

File Descriptions:

1. julia/ARM.jl: The full implementation of the model
2. julia/armcalibrate.jl: Implementations of approximate bayesian calibration and the adoptive metropolis hastings sampling
3. julia/armoptimize.jl: Support utilities for fitting the model to ANES data and displaying the results
4. julia/exploreParams.jl: Performs a parameter grid search to better understand the models sensitivity
5. julia/fitARM.jl: Finds the optimal parameters for each of the four versions of the model, relative to ANES data
6. julia/ARM-fit.ipynb: Explores the model variants' dynamics relative to ANES data
7. julia/ARM-calibrate.ipynb: Tests bayesian parameter estimation methods
8. julia/ARM-develop.ipynb: Explores the dynamics of the original model across different parameterizations

References:

Axelrod, Robert, Joshua J. Daymude, and Stephanie Forrest. "Preventing extreme polarization of political attitudes." Proceedings of the National Academy of Sciences 118.50 (2021)
