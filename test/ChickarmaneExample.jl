# Chickarmane example
using Optim
using DifferentialEquations
using Plots
gr()

include("MinimumActionPath.jl");    using MAP
include("ChickarmaneSystem.jl")

# Specify the parameters
Tspan = 500.0;
N = 150;

xₚ = [40.0 80.0 80.0 80.0]; # Artificial state on the high "plateau"

# Use Optim to optimise a single path φ
resObj = MAP.MAP_Opt(f, g, xₚ, xₑ, Tspan,N, res);
res = Optim.minimizer(resObj)
plot(res[:,4], res[:,1], xlims=[0,100],ylims=[0,80])
