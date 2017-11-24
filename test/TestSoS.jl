# Run a series of tests of the NormDecomp function

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials
# using SCS
# using Mosek
using CSDP

include("../src/NormalSoS.jl")
using NormalSoS

#= Notes on solvers
SCS    - Tends to be very slow, sometimes gives optimal solutions
Mosek  - Cannot run on some problems (ranged constraints)
CSDP   - Faster than SCS but requires more iterations or higher order bounding
         polynomial to achieve optimality
=#

## Example 1: Quadratic system from Zhou et al (2012) - Y
@polyvar x[1:2]
F1(x::Vector) = [-x[1] + 2.0x[2]^2;
     -x[1]*x[2] - 2.0x[2]];
f1 = F1(x);
@time Ueg1 = NormalSoS.normdecomp(f1,x, CSDPSolver())
# PlotLandscape(f,Ueg1,vars,[-3 3],[-3 3]);
# PlotVectors(f,Ueg1,vars,[-3 3],[-3 3]);
# CheckNorm(f,Ueg1,vars)


## Example 2: Quartic system from Zhou et al (2012) - Y
@polyvar x[1:2]
F2 = [-1.0 + 9.0x[1] - 2.0x[1]^3 + 9.0x[2] - 2.0x[2]^3;
      1.0 - 11.0x[1] + 2.0x[1]^3 + 11.0x[2] - 2.0x[2]^3];
f2 = F2(x);
@time Ueg2 = NormalSoS.normdecomp(f2,x, CSDPSolver(),1,4)


## Example 3: One-dimensional bistable system - Y

# @polyvar x[1]
# f3 = [x[1] - x[1]^3 + 0.1];
# @time Ueg3 = normalSoS.NormDecomp(f3,x)


## Example 4: Two-dimensional bistable system with curl dynamics - Y
α = 0.5;    λ = 0.5α;    β = -0.05;    c = 0.0;
@polyvar x[1:2]
F4 = [2α*x[1] - 4λ*x[1]^3 - β + 4c*λ*x[2]^3;
     2c*α*x[1] - 4c*λ*x[1]^3 - c*β - 4λ*x[2]^3];
f4 = F4(x);
@time Ueg4 = NormalSoS.normdecomp(f4,x, CSDPSolver(),1,4)


## Example 6: Fei's 4dof Michaelis-Menten enzyme dynamics model
# Unable to find something reasonable so converges to very small
# coefficients. This is also true of the standard Lyapunov method.

@polyvar x[1:4]
h1 = 2.0;    h2 = 0.2;    h3 = 3.0;
F6 = [-h1*x[1]x[2] + h2*x[3];
     -h1*x[1]x[2] + (h2+h3)x[3];
     h1*x[1]x[2] - (h2+h3)x[3];
     h3*x[3]];
f6 = F6(x);
Ueg6 = NormalSoS.normdecomp(f6,x, CSDPSolver())
