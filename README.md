# SDEtools

This repository provides a set of tools for the analysis of stochastic differential equations, with particular application to biological systems.

## Potential landscape generation
It is popular for an SDE to be expressed in terms of a potential function, such that the deterministic part is expressed in terms of the gradient of this potential. In general, this involves decomposing the deterministic part into a gradient and curl contribution. In this toolbox we attempt to perform this decomposition such that these components are orthogonal to each other. The method is based upon the sum of squares (SoS) approach to generating a Lyapunov function, and is impemented in both Julia and Matlab.

## Minimum action paths
Related to the concept of potential landscapes is the minimum action path (MAP), that determines the most probable path between any two points in the state space. For a pure gradient system this is determined by the potential function, while for a combined system the path depends on both components. We compute the MAP using the nudged elastic band (NEB) method.
