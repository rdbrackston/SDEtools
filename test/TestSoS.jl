# Run a series of tests of the NormDecomp function

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials
using SCS, Mosek, CSDP
using Plots
gr()

include("../src/NormalSoS.jl")
using NormalSoS

#= Notes on solvers
SCS    - Tends to be very slow, sometimes gives optimal solutions.
Mosek  - Fastest. Formerly could not run on some problems (ranged constraints).
CSDP   - Faster than SCS but requires more iterations or higher order bounding
         polynomial to achieve optimality.
=#

## Example 1: Quadratic system from Zhou et al (2012) - Y
@polyvar x1[1:2]
F1(x::Vector) = [-x[1] + 2.0x[2]^2;
     -x[1]*x[2] - 2.0x[2]];
f1 = F1(x1);
@time Ueg1 = NormalSoS.normdecomp(f1,x1, MosekSolver())
plt1 = NormalSoS.plotlandscape(f1,Ueg1,x1,([-3 3],[-3 3]),false);    plot(plt1)
NormalSoS.checknorm(f1,Ueg1,x1)
plttmp = NormalSoS.plotvectors(f1,x1,([-3 3],[-3 3]));    plot(plttmp)


## Example 2: Quartic system from Zhou et al (2012) - Y
@polyvar x2[1:2]
F2(x::Vector) = [-1.0 + 9.0x[1] - 2.0x[1]^3 + 9.0x[2] - 2.0x[2]^3;
      1.0 - 11.0x[1] + 2.0x[1]^3 + 11.0x[2] - 2.0x[2]^3];
f2 = F2(x2);
@time Ueg2 = NormalSoS.normdecomp(f2,x2, MosekSolver(),1,4)
plt2 = NormalSoS.plotlandscape(f2,Ueg2,x2,([-3 3],[-3 3]));    plot(plt2)
NormalSoS.checknorm(f2,Ueg2,x2)


## Example 3: One-dimensional bistable system - Y
@polyvar x3;    x3 = [x3];
F3(x::Vector) = [x[1] - x[1]^3 + 0.1];
f3 = F3(x3);
@time Ueg3 = NormalSoS.normdecomp(f3,x3, MosekSolver(),1,4)
NormalSoS.checknorm(f3,Ueg3,x3)


## Example 4: Two-dimensional bistable system with curl dynamics - Y
α = 0.5;    λ = 0.5α;    β = -0.05;    c = 0.5;
@polyvar x4[1:2]
F4(x::Vector) = [2α*x[1] - 4λ*x[1]^3 - β + 4c*λ*x[2]^3;
     2c*α*x[1] - 4c*λ*x[1]^3 - c*β - 3λ*x[2]^3];
f4 = F4(x4);
@time Ueg4 = NormalSoS.normdecomp(f4,x4, MosekSolver(),1,4)
plt4 = NormalSoS.plotlandscape(f4,Ueg4,x4,([-3 3],[-3 3]),true);    plot(plt4)
NormalSoS.checknorm(f4,Ueg4,x4)


## Example 5:
@polyvar x5[1:2]
F5(x::Vector) = [-x[1] + x[2] + x[1]x[2];
                 -x[1] - x[1]^2];
f5 = F5(x5);
@time Ueg5 = NormalSoS.normdecomp(f5,x5, MosekSolver(),1,4)
plt4 = NormalSoS.plotlandscape(f4,Ueg4,x4,([-3 3],[-3 3]));    plot(plt4)
NormalSoS.checknorm(f4,Ueg4,x4)


## Example 7: From Papachristodolou and Prajna (2005)
# Finds a reasonable landscape and gets close to orthogonality
@polyvar x7[1:4]
F7(x::Vector) = [-x[1] + x[2]^3 - 3*x[3]*x[4];
                 -x[1] - x[2]^3;
                 x[1]*x[4] - x[3];
                 x[1]*x[3] - x[4]^3];
f7 = F7(x7);
@time Ueg7 = NormalSoS.normdecomp(f7,x7, MosekSolver(),1,4)
plt7 = NormalSoS.plotlandscape(f7,Ueg7,x7,([-3 3],[-3 3]), true);    plot(plt7)
NormalSoS.checknorm(f7,Ueg7,x7)


##
