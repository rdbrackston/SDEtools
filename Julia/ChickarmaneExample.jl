# Chickarmane example
using Optim
using DifferentialEquations
using Gadfly

include("MinimumActionPath.jl")
using MAP

# Specify the equations of the model
(k0, c0, c1, c2, c3, c4, e0, e1, e2, a0, a1, a2, b0, b1, b2, b3, δ) = [0.005, 0.01, 0.4, 1.0, 0.1, 0.00135, 0.01, 1.0, 1.0, 0.01, 1.0, 5.0, 0.005, 0.005, 1.0, 1.0, 0.01];
α, LIF, I3 = [0.0, 0.0, 6.0];

r1(t,u) = k0*u[2]*(c0 + c1*u[1]*u[1] + k0*u[2] + c2*LIF) /
         ( 1 + k0*u[2]*(c1*u[1]*u[1] + k0*u[2] + c2*LIF + c3*u[3]*u[3]) +
                       + c4*u[2]*u[4]*u[4] );
r2(t,u) = δ*u[1];
r3(t,u) = α + (e0 + e1*u[2]) / (1 + e1*u[2] + e2*u[4]*u[4]);
r4(t,u) = δ*u[2];
r5(t,u) = (a0 + a1*u[2]) / (1 + a1*u[2] + a2*I3);
r6(t,u) = δ*u[3];
r7(t,u) = (b0 + b1*u[4]*u[4] + b3*u[2]) /
          (1 + b1*u[4]*u[4] + b2*u[1]*u[1] + b3*u[2]);
r8(t,u) = δ*u[4];

f(X::Vector) = [r1(0,X) - r2(0,X),
                r3(0,X) - r4(0,X),
                r5(0,X) - r6(0,X),
                r7(0,X) - r8(0,X)];
g(X::Vector) = [sqrt(r1(0,X) + r2(0,X)),
                sqrt(r3(0,X) + r4(0,X)),
                sqrt(r5(0,X) + r6(0,X)),
                sqrt(r7(0,X) + r8(0,X))];

cond_ss = [83.0,97.0,76.0,1.0]; # Stem cell state
I3 = 2.0;
ss_problem = DifferentialEquations.SteadyStateProblem(
                    (t, X) -> [r1(t, X) - r2(t, X),
                              r3(t, X) - r4(t, X),
                              r5(t, X) - r6(t, X),
                              r7(t, X) - r8(t, X)],
                              cond_ss)
sol = DifferentialEquations.solve(ss_problem);
