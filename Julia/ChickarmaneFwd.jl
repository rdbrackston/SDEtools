# Evaluate optimal forward paths
using Optim
using DifferentialEquations
using Plots
gr()

include("MinimumActionPath.jl");    using MAP
include("ChickarmaneSystem.jl")

# Specify the parameters
N = 150;

# Perform the optimisation
Tvec = [];    Svec = [];
resObj = MAP.T_Opt!(Tvec,Svec, f,g,xₛ,xₑ,(500.,5000.),N)
res = Optim.minimizer(resObj)
plt = plot(Tvec,Svec)

resFull = hcat(collect(linspace(0,Tvec[end],N-2)),resrev);
writedlm("ChickFwdPth.txt",resFull);
