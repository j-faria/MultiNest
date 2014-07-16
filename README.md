**IMPORTANT** MultiNest is not my (João Faria) project; all copyright belongs to its rightful owners (see below).


MultiNest
=========

Efficient and Robust Bayesian Inference
---------------------------------------


MultiNest is a Bayesian inference tool which calculates the evidence and explores the parameter space which may 
contain multiple posterior modes and pronounced (curving) degeneracies in moderately high dimensions.
 
 
MultiNest is freely available for academic use, subject to some [licence](./MultiNest_License.pdf) restrictions. 
Please contact [Farhan Feroz](http://www.mrao.cam.ac.uk/~ff235/) for commercial licence.

 
MultiNest can be downloaded from the files section of [this website](http://ccpforge.cse.rl.ac.uk/gf/project/multinest/).
 
 
 
If you use MultiNest, please acknowledge the following MultiNest papers:

[arXiv:0809.3437](http://xxx.lanl.gov/abs/0809.3437)
[arXiv:0704.3704](http://xxx.lanl.gov/abs/0704.3704)
[arXiv:1306.2144](http://xxx.lanl.gov/abs/1306.2144)
 
 

Language Support
----------------

MultiNest is written in Fortran 90 but wrappers & in the case of Matlab a native language implementation, 
is available for likelihood and prior functions written in following languages:

C/C++: wrappers provided with MultiNest package
R: RMultiNest, wrapper for MultiNest written by Johannes Buchner.
Python: PyMultiNest, wrapper for MultiNest written by Johannes Buchner.
Matlab: MatlabMultiNest written by Matthew Pitkin & Joe Romano, is a Matlab implementation of the algorithm 
        described in arXiv:0809.3437, although it lacks a few features included in the original Fortran 
        version of MultiNest. It can be downloaded from the files section of 
        [this website](http://ccpforge.cse.rl.ac.uk/gf/project/multinest/).
Read the MultiNest README file for more details.
 

Plotting & Visualisation
------------------------

GetDist: Output files produced by MultiNest are compatible with GetDist? plotting package which is part of CosmoMC.
SuperPlot: written by Andrew Fowlie, can calculate and plot posterior distributions, profile likelihoods, 
           confidence intervals, credible regions etc. from MultiNest output files.
PyMultiNest: has support for analysing the output from MultiNest & monitoring the progress.
 

Packages Using MultiNest
------------------------

BAMBI: Blind Accelerated Multimodal Bayesian Inference
Bayes-X is a Bayesian inference tool for the analysis of X-ray observations of clusters of galaxies
BXA: Bayesian X-ray analysis (nested sampling for Xspec and Sherpa)
Cosmo++: An Object Oriented C++ Library for Cosmology
ModeCode: Bayesian Parameter Estimation and Model Comparison for Inflation
Monte Python: Monte Carlo Code for Cosmic Linear Anisotropy Solving System (CLASS) in Python
pSNid III is a software package to probabilistically identify/classify supernova light curves
SIMTOI: SImulation and Modeling Tool for Optical Interferometry
SuperBayeS for fast and efficient analysis of supersymmetry theories in particle physics.
SuperPy: Scan the CMSSM with (Py)MultiNest and plot the results with a new plotting GUI.
TempoNest is a Bayesian software package for the analysis of pulsar timing data.
