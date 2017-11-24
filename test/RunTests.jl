# Run a series of tests of the NormDecomp function

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials
# using SCS
using Mosek
include("normal-SoS.jl")
using NormalSoS


## Example 1: Quadratic system from Zhou et al (2012) - Y

@polyvar x[1:2]
f1 = [-x[1] + 2.0x[2]^2;
     -x[1]*x[2] - 2.0x[2]];
@time Ueg1 = normalSoS.normdecomp(f1,x, MosekSolver())
# PlotLandscape(f,Ueg1,vars,[-3 3],[-3 3]);
# PlotVectors(f,Ueg1,vars,[-3 3],[-3 3]);
# CheckNorm(f,Ueg1,vars)


## Example 2: Quartic system from Zhou et al (2012) - Y

@polyvar x[1:2]
f2 = [-1.0 + 9.0x[1] - 2.0x[1]^3 + 9.0x[2] - 2.0x[2]^3;
      1.0 - 11.0x[1] + 2.0x[1]^3 + 11.0x[2] - 2.0x[2]^3];
@time Ueg2 = normalSoS.normdecomp(f2,x)


## Example 3: One-dimensional bistable system - Y

# @polyvar x[1]
# f3 = [x[1] - x[1]^3 + 0.1];
# @time Ueg3 = normalSoS.NormDecomp(f3,x)


## Example 4: Two-dimensional bistable system with curl dynamics
# Works well, only if the curl component is significant (~0.1), otherwise
# it generates something non orthogonal. Using a 4th-order lower bound
# polynomial helps.

@polyvar x[1:2]
f4 = [-x[1]^3 + x[1] - 2.0x[1]^3*x[2]^4 + 2.0x[1]*x[2]^4;
     -4.0x[2]^3 - 2.0x[1]^4*x[2]^3 + 4.0x[1]^2*x[2]^3];
fc4 = [4.0x[2]^3 + 2.0x[1]^4*x[2]^3 - 4.0x[1]^2*x[2]^3;
         x[1] - x[1]^3 - 2.0x[1]^3*x[2]^4 + 2.0x[1]*x[2]^4];
#f4 = f4 + 0.1fc4;
Ueg4 = normalSoS.normdecomp(f4,x)


## Example 6: Fei's 4dof Michaelis-Menten enzyme dynamics model
# Unable to find something reasonable so converges to very small
# coefficients. This is also true of the standard Lyapunov method.

@polyvar x[1:4]
h1 = 2.0;    h2 = 0.2;    h3 = 3.0;
f6 = [-h1*x[1]x[2] + h2*x[3];
     -h1*x[1]x[2] + (h2+h3)x[3];
     h1*x[1]x[2] - (h2+h3)x[3];
     h3*x[3]];

Ueg6 = normalSoS.normdecomp(f6,x)
