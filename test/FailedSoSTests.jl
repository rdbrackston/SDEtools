# Script containing failing examples

## Example 5: Fei's 4dof Michaelis-Menten enzyme dynamics model
# Unable to find something reasonable. Issue may be because x[3] does not
# feature in the state equations.
@polyvar x5[1:4]
h1 = 0.2;    h2 = 0.02;    h3 = 0.3;
F5(x::Vector) = [-h1*x[1]x[2] + h2*x[3];
     -h1*x[1]x[2] + (h2+h3)x[3];
     h1*x[1]x[2] - (h2+h3)x[3];
     h3*x[3]];
f5 = F5(x5);
Ueg5 = NormalSoS.lyapunov(f5,x5,MosekSolver(),2,true)
basis = NormalSoS.minimalbasis(f5,x5);    Ueg5 = NormalSoS.normopt2(f5,x5,basis,MosekSolver(),4)
@time Ueg5 = NormalSoS.normdecomp(f5,x5, MosekSolver(),1,4)
plt5 = NormalSoS.plotlandscape(f5,Ueg5,x5,([0 3],[0 3]));    plot(plt5)
NormalSoS.checknorm(f5,Ueg5,x5)


## Example 5: More comlicated two-dimensional bistable system with curl dynamics
α = 0.5;    λ = 0.5α;    β = -0.05;    c = 0.5;    d = 1.0;
@polyvar x5[1:2]
F5(x::Vector) = [-4d*λ*x[1]^3 + 2d*α*x[1] - d*β - 4λ^2*x[1]^3*x[2]^4 + 2α*λ*x[1]*x[2]^4 - β*λ*x[2]^4;
                 -4d*λ*x[2]^3 - 4λ^2*x[1]^4*x[2]^3 + 4α*x[1]^2*x[2]^3 - 4β*λ*x[1]*x[2]^3];
F5c(x::Vector) = [4d*λ*x[2]^3 + 4λ^2*x[1]^4*x[2]^3 - 4α*x[1]^2*x[2]^3 + 4β*λ*x[1]*x[2]^3;
                 -4d*λ*x[1]^3 + 2d*α*x[1] - d*β - 4λ^2*x[1]^3*x[2]^4 + 2α*λ*x[1]*x[2]^4 - β*λ*x[2]^4];
f5 = F5(x5) + c*F5c(x5);
Uan5 = (d + λ*x5[1]^4 - α*x5[1]^2 + β*x5[1])*(d + λ*x5[2]^4);
basis = NormalSoS.minimalbasis(f5,x5);    Ueg5 = NormalSoS.normopt2(f5,x5,basis,MosekSolver(),4)
@time Ueg5 = NormalSoS.normdecomp(f5,x5, MosekSolver(),1,4,:minimal)#,Uan5)
plt5 = NormalSoS.plotlandscape(f5,Ueg5,x5,([-3 3],[-3 3]),false);    plot(plt5)
NormalSoS.checknorm(f5,Ueg5,x5)


## Example 6: Fei's 10dof Michaelis-Menten enzyme dynamics model
# Failed - memory requirements too large.
@polyvar x6[1:10]
(h1, h2,  h3, h4,  h5,   h6,  h7, h8,    h9, h10, h11,    h12, h13,   h14, h15,  h16,   h17,    h18,  h19,  h20, h21,  h22,  h23,   h24) =
(30.,6e-5,30.,6e-5,1.221,6e-5,5.4,0.0048,30.,6e-5,9.24e-5,0.99,0.0168,1.35,0.075,0.2448,0.00678,0.018,0.012,11.1,0.075,0.828,0.0072,0.2442);
F6(x::Vector) = [-(h17+h18)x[1] + h2*x[3] + h15*x[4] + h16*x[10] - h1*x[1]x[2] - h14*x[1]x[6];
                 -h7*x[2] + (h2+h6)x[3] + (h4+h5)x[5] + h8*x[7] - h1*x[1]x[2] - h3*x[2]x[4];
                 -(h2+h6)x[3] + h21*x[5] + h22*x[9] + h1*x[1]x[2] - h20*x[3]x[6];
                 -(h15+h24)x[4] + h4*x[5] + h14*x[1]x[6] - h3*x[2]x[4];
                 -(h4+h5+h21)x[5] + h3*x[2]x[4] + h20*x[3]x[6];
                 (h15+h24)x[4] + (h5+h21)x[5] - h23*x[6] - h14*x[1]x[6] - h20*x[3]x[6];
                 h7*x[2] - h8*x[7] + h10*x[9] - h9*x[7]x[8];
                 h18*x[1] - h19*x[8] + h10*x[9] - h9*x[7]x[8];
                 -(h10+h22)x[9] + h9*x[7]x[8];
                 h11 - h13*x[10] + h12*x[7]^2 ];
f6 = F6(x6);
@time Ueg6 = NormalSoS.normdecomp(f6,x6, CSDPSolver(),1,4)
plt6 = NormalSoS.plotlandscape(f6,Ueg6,x6,([-3 3],[-3 3]));    plot(plt6)
NormalSoS.checknorm(f6,Ueg6,x6)


## Example 8: From Papachristodolou and Prajna (2003)
# Failed - memory requirements too large.
@polyvar x8[1:6]
F8(x::Vector) = [-x[1]^3 + 4*x[2]^3 - 6*x[3]*x[4];
                 -x[1] - x[2] + x[5]^3;
                 x[1]x[4] - x[3] + x[4]x[6];
                 x[1]x[3] + x[3]x[6] - x[4]^3;
                 -2*x[2]^3 - x[5] + x[6];
                 -3x[3]x[4] - x[5]^3 - x[6]];
f8 = F8(x8);
@time Ueg8 = NormalSoS.normdecomp(f8,x8, CSDPSolver(),1,4)
plt8 = NormalSoS.plotlandscape(f8,Ueg8,x8,([-3 3],[-3 3]));    plot(plt8)
NormalSoS.checknorm(f8,Ueg8,x8)


## 8: From Papachristodolou and Prajna (2005)
# Fails to find anything sensible - possibly due to limit cycle or unstable behaviour.
a = 0.1;    b = 0.4;
@polyvar x8[1:2]
F8(x::Vector) = [a - x[1] + x[1]^2*x[2];
                 b - x[1]^2*x[2]];
f8 = F8(x8);
@time Ueg8 = NormalSoS.normdecomp(f8,x8, CSDPSolver(),1,2)
plt6 = NormalSoS.plotlandscape(f6,Ueg6,x6,([-2 2],[-2 2]),true);    plot(plt6)
NormalSoS.checknorm(f6,Ueg6,x6)


# Old positive definiteness constraints
# @polyconstraint m V ≥ ϵ[1]*x[1]^o+ϵ[2]*x[2]^o;
# bnd = monomials(x,collect(2:2:o), m -> exponents(m)[1]!=1 && exponents(m)[1]!=3);
# bnd = monomials(x,collect(2:2:o), m -> exponents(m)[1]==0 || exponents(m)[2]==0)
# bnd = monomials(x,[o], m -> exponents(m)[1]!=1 && exponents(m)[1]!=3);
# bnd = monomials(x,[o], m -> exponents(m)[1]==0 || exponents(m)[2]==0);
