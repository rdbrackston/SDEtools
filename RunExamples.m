%% Script to run through each of the example systems and plot the outputs

clearvars; clc; echo off;


%% Example 1: Quadratic system from Zhou et al (2012) - Y

syms x1 x2
vars = [x1;x2];

f = [-1.0*x1 + 2*x2^2;
     -x1*x2 - 2*x2];

% Ueg1 = NormDecomp(f,vars)
Ueg1 = Lyapunov(f,vars,2)
% PlotLandscape(f,Ueg1,vars,[-3 3],[-3 3]);
% PlotVectors(f,Ueg1,vars,[-3 3],[-3 3]);
CheckNorm(f,Ueg1,vars)


%% Example 2: Quartic system from Zhou et al (2012) - Y

syms x1 x2;
vars = [x1; x2];

% Constructing the vector field dx/dt = f
f = [-1 + 9*x1 - 2*x1^3 + 9*x2 - 2*x2^3;
     1 - 11*x1 + 2*x1^3 + 11*x2 - 2*x2^3];

Ueg2 = NormDecomp(f,vars)
% PlotLandscape(f,Ueg2,vars,[-3 3],[-3 3])
% PlotVectors(f,Ueg2,vars,[-3 3],[-3 3])
CheckNorm(f,Ueg2,vars)


%% Example 3: One-dimensional bistable system - Y

syms x1
vars = x1;

f = x1 - x1^3 + 0.1;

Ueg3 = NormDecomp(f,vars)
% PlotLandscape(f,Ueg3,vars,[-2 2],[-2 2]);
CheckNorm(f,Ueg3,vars)


%% Example 4: Two-dimensional bistable system with curl dynamics
% Works well, only if the curl component is significant (~0.1), otherwise
% it generates something non orthogonal. Using a 4th-order lower bound
% polynomial helps.

syms x1 x2
vars = [x1;x2];

f = [-x1^3 + x1 - 2*x1^3*x2^4 + 2*x1*x2^4;
     -4*x2^3 - 2*x1^4*x2^3 + 4*x1^2*x2^3];
fPerp = [4*x2^3 + 2*x1^4*x2^3 - 4*x1^2*x2^3;
         x1 - x1^3 - 2*x1^3*x2^4 + 2*x1*x2^4];
f = f + 0.1*fPerp;

Ueg4 = NormDecomp(f,vars,2,4)
% PlotLandscape(f,Ueg4,vars,[-1.5 1.5],[-1.5 1.5]);
% PlotVectors(f,Ueg4,vars,[-2 2],[-1.5 1.5]);
CheckNorm(f,Ueg4,vars)


%% Example 5 Lotka-Volterra
% Unable to find something reasonable so converges to very small
% coefficients. This is also true of the standard Lyapunov method.

syms x1 x2
vars = [x1; x2];

a=2/3;    b=4/3;    c=1;    d=2;
f = [a*x1 - b*x1*x2;
     c*x1*x2 - d*x2];

% Ueg5 = NormDecomp(f,vars,2)
Ueg5 = Lyapunov(f,vars,2)
% PlotLandscape(f,Ueg5,vars,[0 20],[0 20]);
% PlotVectors(f,Ueg5,vars,[0 20],[0 20]);
CheckNorm(f,Ueg5,vars)


%% Example 6: Fei's 4dof Michaelis-Menten enzyme dynamics model
% Unable to find something reasonable so converges to very small
% coefficients. This is also true of the standard Lyapunov method.

syms x1 x2 x3 x4
vars = [x1;x2;x3;x4];

h1 = 2;    h2 = 0.2;    h3 = 3;
f = [-h1*x1*x2 + h2*x3;
     -h1*x1*x2 + (h2+h3)*x3;
     h1*x1*x2 - (h2+h3)*x3;
     h3*x3];

% Ueg6 = NormDecomp(f,vars)
Ueg6 = Lyapunov(f,vars,2)
% PlotLandscape(f,Ueg6,vars,[0 10],[0 10]);
% PlotVectors(f,Ueg6,vars,[0 10],[0 10]);
CheckNorm(f,Ueg6,vars)


%% Example 7: From Papachristodolou and Prajna (2005)
% Finds a reasonable loking Lyapunov function but the decomposition is not
% orthogonal

syms x1 x2 x3 x4
vars = [x1;x2;x3;x4];

f = [-x1 + x2^3 - 3*x3*x4;
     -x1 - x2^3;
      x1*x4 - x3;
      x1*x3 - x4^3];

Ueg7 = NormDecomp(f,vars)
PlotLandscape(f,Ueg7,vars,[-3 3],[-3 3]);
PlotVectors(f,Ueg7,vars,[-3 3],[-3 3]);
CheckNorm(f,Ueg7,vars)


