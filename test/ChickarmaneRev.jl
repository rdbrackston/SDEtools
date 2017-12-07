# Evaluate optimal reverse paths
using Optim
using DifferentialEquations
using Plots
gr()

include("../src/MinimumActionPath.jl");    using MAP
include("ChickarmaneSystem.jl")

# Specify the parameters
N = 150;


################################################################################
# Perform the full temporal optimisation
trackerObj = MAP.OptimisationTracker([],[],[],[]);
resObj = MAP.optimaltime!(trackerObj, f,g,xₑ,xₛ,(1000.,5000.),N);
plt = plot(trackerObj.pointVec[4:end],trackerObj.valueVec[4:end], line=0,marker=2);
T = trackerObj.pointVec[end];


################################################################################
# Alternatively perform a fixed time optimisation
T = 2500;   pth = MAP.makepath(xₑ,xₛ,N);
conv = false;
while !conv
    resObj = MAP.optimalpath(f,g,xₑ,xₛ,T,N, pth);
    pth = Optim.minimizer(resObj);
    conv = Optim.converged(resObj);
end


################################################################################
# Plot path and save to file
res = vcat(xₑ,Optim.minimizer(resObj),xₛ);
plt = plot(res[:,4],res[:,1], line=0,marker=4)

resFull = hcat(collect(linspace(0,T,N)),res);
plot!(resFull[:,1],resFull[:,2], label="Nanog")
plot!(resFull[:,1],resFull[:,3], label="Oct4")
plot!(resFull[:,1],resFull[:,4], label="Fgf4")
plot!(resFull[:,1],resFull[:,5], label="Gata6")

writedlm("ChickRevPth2.txt",resFull);
