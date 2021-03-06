# Run a series of tests of the NormDecomp function

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials
using Mosek
using Plots
gr()

include("../src/NormalSoS.jl")
using NormalSoS


## Example 1: Quadratic system from Zhou et al (2012)
@polyvar x1[1:2]
F1(x::Vector) = [-x[1] + 1.0x[2]^2;
     -1.5*x[1]*x[2] - 2.0x[2]];
f1 = F1(x1);
Ueg1 = NormalSoS.normdecomp(f1,x1, MosekSolver(),2, :minimal)
plt1 = NormalSoS.plotlandscape(f1,Ueg1,x1,([-3 3],[-3 3]),true);    plot(plt1)
NormalSoS.checknorm(f1,Ueg1,x1)


## Example 2: Quartic system from Zhou et al (2012)
@polyvar x2[1:2]
F2(x::Vector) = [-1.0 + 9.0x[1] - 2.0x[1]^3 + 9.0x[2] - 2.0x[2]^3;
      1.0 - 11.0x[1] + 2.0x[1]^3 + 11.0x[2] - 2.0x[2]^3];
f2 = F2(x2);
@time Ueg2 = NormalSoS.normdecomp(f2,x2, MosekSolver(),1, :minimal)
plt2 = NormalSoS.plotlandscape(f2,Ueg2,x2,([-3 3],[-3 3]),true);    plot(plt2)
NormalSoS.checknorm(f2,Ueg2,x2)


## Example 3: One-dimensional bistable system
@polyvar x3;    x3 = [x3];
F3(x::Vector) = [x[1] - x[1]^3 + 0.1];
f3 = F3(x3);
@time Ueg3 = NormalSoS.normdecomp(f3,x3, MosekSolver(),1, :minimal)
NormalSoS.checknorm(f3,Ueg3,x3)


## Example 4: Two-dimensional bistable system with curl dynamics
α = 0.5;    λ = 0.5α;    β = -0.05;    c = 0.5;
@polyvar x4[1:2]
F4(x::Vector) = [2α*x[1] - 4λ*x[1]^3 - β + 4c*λ*x[2]^3;
     2c*α*x[1] - 4c*λ*x[1]^3 - c*β - 4λ*x[2]^3];
f4 = F4(x4);
Uan4 = λ*(x4[1]^4+x4[2]^4) - α*x4[1]^2 + β*x4[1];
@time Ueg4 = NormalSoS.normdecomp(f4,x4, MosekSolver(),1, :minimal)
plt4 = NormalSoS.plotlandscape(f4,Ueg4,x4,([-3 3],[-3 3]),true);    plot(plt4)
NormalSoS.checknorm(f4,Ueg4,x4)


## Example 5: Three-dimensional bistable system
α = 0.5;    λ = 0.5α;    β = -0.05;    c = 0.0;
@polyvar x5[1:3]
F5(x::Vector) = [2α*x[1] - 4λ*x[1]^3 - β;
                 -4λ*x[2]^3;
                 -4λ*x[3]^3];
f5 = F5(x5);
@time Ueg5 = NormalSoS.normdecomp(f5,x5, MosekSolver(),0, :minimal)
plt5 = NormalSoS.plotlandscape(f5,Ueg5,x5,([-3 3],[-3 3]),false);    plot(plt5)
NormalSoS.checknorm(f5,Ueg5,x5)


## Example 6: The Maier-Stein Model
γ = 1.0;    μ = 1.0;
@polyvar x6[1:2]
F6(x::Vector) = [x[1] - x[1]^3 - γ*x[1]x[2]^2;
                 -μ*(x[1]^2 + 1)x[2]];
f6 = F6(x6);
Uan6 = -0.5*x6[1]^2 + 0.25*x6[1]^4 + 0.5γ*x6[2]^2 + 0.5γ*x6[1]^2*x6[2]^2;
@time Ueg6 = NormalSoS.normdecomp(f6,x6, MosekSolver(),2, :minimal)
plt6 = NormalSoS.plotlandscape(f6,Ueg6,x6,([-2 2],[-2 2]),true);    plot(plt6)
NormalSoS.checknorm(f6,Ueg6,x6)


## Example 7: From Papachristodolou and Prajna (2005)
@polyvar x7[1:4]
F7(x::Vector) = [-x[1] + x[2]^3 - 3*x[3]*x[4];
                 -x[1] - x[2]^3;
                 x[1]*x[4] - x[3];
                 x[1]*x[3] - x[4]^3];
f7 = F7(x7);
@time Ueg7 = NormalSoS.normdecomp(f7,x7, MosekSolver(),0,:minimal)
plt7 = NormalSoS.plotlandscape(f7,Ueg7,x7,([-3 3],[-3 3]), false);    plot(plt7)
NormalSoS.checknorm(f7,Ueg7,x7)


## Example 8: Nonlinear system where extended basis performs better
@polyvar x8[1:2]
F8(x::Vector) = [-x[1]^3 - x[1]*x[2]^2 - x[2]^3 - x[1];
                 -x[1]*x[2]^2 - x[2]^3];
f8 = F8(x8);
Ueg8 = NormalSoS.normdecomp(f8,x8, MosekSolver(),2,:minimal)
NormalSoS.checknorm(f8,Ueg8,x8)
Ueg8 = NormalSoS.normdecomp(f8,x8, MosekSolver(),2,:extended)
NormalSoS.checknorm(f8,Ueg8,x8)


## Example 9: Separable system from Cameron (2012)
@polyvar x9[1:2]
F9(x::Vector) = [1 - x[1] - x[2]^2;
                 x[2] - x[2]^3];
f9 = F9(x9);
@time Ueg9 = NormalSoS.normdecomp(f9,x9, MosekSolver(),1, :minimal)
plt9 = NormalSoS.plotlandscape(f9,Ueg9,x9,([-3 3],[-3 3]), true);    plot(plt9)
NormalSoS.checknorm(f9,Ueg9,x9)
