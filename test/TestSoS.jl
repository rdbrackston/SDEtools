# Run a series of tests of the NormDecomp function

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials
using SCS
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
NormalSoS.plotlandscape(f1,Ueg1,x,([-3 3],[-3 3]));
NormalSoS.checknorm(f1,Ueg1,x)


## Example 2: Quartic system from Zhou et al (2012) - Y
@polyvar x[1:2]
F2(x::Vector) = [-1.0 + 9.0x[1] - 2.0x[1]^3 + 9.0x[2] - 2.0x[2]^3;
      1.0 - 11.0x[1] + 2.0x[1]^3 + 11.0x[2] - 2.0x[2]^3];
f2 = F2(x);
@time Ueg2 = NormalSoS.normdecomp(f2,x, CSDPSolver(),1,4)
plt = NormalSoS.plotlandscape(f2,Ueg2,x,([-3 3],[-3 3]));    plot(plt)
NormalSoS.checknorm(f2,Ueg2,x)


## Example 3: One-dimensional bistable system - Y
@polyvar x;    x = [x];
F3(x::Vector) = [x[1] - x[1]^3 + 0.1];
f3 = F3(x);
@time Ueg3 = NormalSoS.normdecomp(f3,x, CSDPSolver(),1,4)
NormalSoS.checknorm(f3,Ueg3,x)


## Example 4: Two-dimensional bistable system with curl dynamics - Y
α = 0.5;    λ = 0.5α;    β = -0.05;    c = 0.0;
@polyvar x[1:2]
F4(x::Vector) = [2α*x[1] - 4λ*x[1]^3 - β + 4c*λ*x[2]^3;
     2c*α*x[1] - 4c*λ*x[1]^3 - c*β - 4λ*x[2]^3];
f4 = F4(x);
@time Ueg4 = NormalSoS.normdecomp(f4,x, CSDPSolver(),1,4)
plt4 = NormalSoS.plotlandscape(f4,Ueg4,x,([-3 3],[-3 3]));    plot(plt4)
NormalSoS.checknorm(f4,Ueg4,x)


## Example 6: Fei's 4dof Michaelis-Menten enzyme dynamics model
# Unable to find something reasonable so converges to very small
# coefficients. This is also true of the standard Lyapunov method.

@polyvar x[1:4]
h1 = 0.2;    h2 = 0.02;    h3 = 0.3;
F6(x::Vector) = [-h1*x[1]x[2] + h2*x[3];
     -h1*x[1]x[2] + (h2+h3)x[3];
     h1*x[1]x[2] - (h2+h3)x[3];
     h3*x[3]];
f6 = F6(x);
@time Ueg6 = NormalSoS.normdecomp(f6,x, CSDPSolver(),1,4)
plt6 = NormalSoS.plotlandscape(f6,Ueg6,x,([-3 3],[-3 3]));    plot(plt6)
NormalSoS.checknorm(f6,Ueg6,x)


## Example 7: From Papachristodolou and Prajna
@polyvar x[1:4]
F7(x::Vector) = [-x[1] + x[2]^3 - 3*x[3]*x[4];
                 -x[1] - x[2]^3;
                 x[1]*x[4] - x[3];
                 x[1]*x[3] - x[4]^3];
f7 = F7(x);
@time Ueg7 = NormalSoS.normdecomp(f7,x, CSDPSolver(),0,4)
