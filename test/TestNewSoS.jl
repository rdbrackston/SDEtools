# Run a series of tests of the NormDecomp function

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials
using Mosek
using Plots
gr()

include("../src/NormalSoS.jl")
using NormalSoS

# A 3D non-normal example - WORKS and agrees with quasipotential
n = 3;    @polyvar x[1:n]
A = [-5.0 0.0 0.2;
     0.0 -1.5 3.0;
     0.5 -5.0 -1.0];
F(X::Vector) = A*X;
f = F(x);
Ul = NormalSoS.normdecomp(f,x, MosekSolver(),2,:extended)
# V = NormalSoS.initialguess(f,x)
Uu = NormalSoS.upperbound(f,x,true,10)


## Example 1: Quadratic system from Zhou et al (2012)
@polyvar x[1:2]
F(x::Vector) = [-1.0 + 9.0x[1] - 2.0x[1]^3 + 9.0x[2] - 2.0x[2]^3;
      1.0 - 11.0x[1] + 2.0x[1]^3 + 11.0x[2] - 2.0x[2]^3];
f = F(x);
Ul = NormalSoS.normdecomp(f,x, MosekSolver(),2, :minimal)
Uu = NormalSoS.upperbound(f,x,false,10)
