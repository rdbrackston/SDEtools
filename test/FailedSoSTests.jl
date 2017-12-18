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
@time Ueg5 = NormalSoS.normdecomp(f5,x5, MosekSolver(),1,4)
plt5 = NormalSoS.plotlandscape(f5,Ueg5,x5,([-3 3],[-3 3]));    plot(plt5)
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
