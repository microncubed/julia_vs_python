# julia_vs_python

## Description
This repo related to the blogpost <https://microncubed.com/benchmarking-matrix-equation-solvers-python-scipy-vs-julia/>. The is to benchmark Python's SciPy methods for solving sparse matrix equations and to compare them to the equivalent methods in the case of the Julia language. The sparse matrix in question describes the 2-d Poission equation when finite differenced. I use the method of manufactured solutions to set the source term for the Poisson equation so that there will be an analytic solution to compare against.

## Files
There are two notebooks 01_poisson_equation_julia.ipynb and 02_poisson_equation_python.ipynb. They are to be run in the suggested numerical order. Each notebook calls functions, from a functions.jl and functions.py respectively. The second notebook performs the data analysis as well as solving the Poisson equation in Python.

The first notebook outputs txt and csv files containing the data from the Julia solutions.

The png files are the only output from the data analysis.