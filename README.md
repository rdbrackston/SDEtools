# SDEtools

This repository provides a set of tools for the analysis of stochastic differential equations. The code currently runs on Julia v0.6.

## Potential landscape generation
It is popular for an SDE to be expressed in terms of a potential function, such that the deterministic part of the SDE depends upon the gradient of this potential. In general, this involves decomposing the deterministic part into a gradient and curl contribution. In this toolbox we attempt to perform this decomposition such that these components are orthogonal to each other. The method is based upon the sum of squares (SoS) approach to generating a Lyapunov function, and involves solving a convex optimization. The details of the method are given in the following paper:

R. D. Brackston, A. Wynn and M. P. H. Stumpf. Construction of quasi-potentials for stochastic dynamical systems: An optimization approach. [arXiv:1805.07273](https://arxiv.org/pdf/1805.07273.pdf).

## Minimum action paths
Related to the concept of potential landscapes is the minimum action path (MAP), that determines the most probable path between any two points in the state space. For a pure gradient system this is determined by the potential function, while for a combined system the path depends on both components. We compute the MAP using a gradient-based optimization algorithm.
