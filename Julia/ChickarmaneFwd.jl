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
resObj = MAP.T_Opt!(Tvec,Svec, f,g,xₛ,xₑ,(100.,4000.),N)
res = vcat(xₛ,Optim.minimizer(resObj),xₑ);

plt = plot(Tvec,Svec, line=0,marker=2)
plt = plot(res[:,4],res[:,1], line=0,marker=4)


resFull = hcat(collect(linspace(0,Tvec[end],N)),res);
writedlm("ChickFwdPth.txt",resFull);
