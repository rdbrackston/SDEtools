# Evaluate optimal reverse paths
using Optim
using DifferentialEquations
using Plots
gr()

include("../src/MinimumActionPath.jl");    using MAP
include("ChickarmaneSystem.jl")

# Specify the parameters
N = 150;

# Perform the optimisation
Tvec = [];    Svec = [];
resObj = MAP.optimaltime!(Tvec,Svec, f,g,xₑ,xₛ,(100.,4000.),N)
res = vcat(xₑ,Optim.minimizer(resObj),xₛ);

plt = plot(Tvec,Svec, line=0,marker=2)
plt = plot(res[:,4],res[:,1], line=0,marker=4)


resFull = hcat(collect(linspace(0,Tvec[end],N)),res);
writedlm("ChickRevPth.txt",resFull);
