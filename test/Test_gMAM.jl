# Test geometric minimum action method

using Optim
using DifferentialEquations
using Plots
gr()

include("../src/GeomMinActPath.jl");    using gMAM
include("../src/MinimumActionPath.jl");    using MAP

α = 0.5;
λ = 0.5α;
β = -0.05;
σ = 0.4;

N = 100;
n=2;
if n==1
    f(x::Vector) = [2α*x[1]-4λ*x[1]^3-β];
    g(x::Vector) = [σ];
    x₀ = [-1.0];  # Start point (row vector)
    xₑ = [1.0];  # End point (row vector)
elseif n==2
    f(x::Vector) = [2α*x[1]-4λ*x[1]^3-β + 4λ*x[2]^3;
                    -4λ*x[2]^3 + 2α*x[1]-4λ*x[1]^3-β];
    g(x::Vector) = σ*[1.0; 1.0];
    x₀ = [-1.0 0.0];  # Start point (row vector)
    xₑ = [1.0 0.0];  # End point (row vector)
else
    print("Invalid n chosen.")
end

φ₀ = [MAP.makepath(x₀,[0.0 1.0],N);MAP.makepath([0.0 1.0],xₑ,N)];
N = length(φ₀[:,1]+2);
S = gMAM.action(φ₀,f,g,x₀,xₑ,N,n)
grad = zeros(φ₀);    gMAM.actiongradient!(grad, φ₀,f,g,x₀,xₑ,N,n)

# Use Optim to optimise over path φ
resObj = gMAM.optimalpath(f, g, x₀, xₑ, N, φ₀, 10);
res = Optim.minimizer(resObj)
plot(res[:,1],res[:,2], marker=2, line=0)
