using Optim
using ForwardDiff
# using ReverseDiff
using Gadfly

include("MinimumActionPath.jl")
using MAP

α = 0.5;
λ = 0.5α;
β = 0.0;
σ = 0.2;

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

Tspan = 18.0;  # Time span
N = 101;
# xᵢ = [0.0, 0.0];  # Intermediate point

# Examine the path
φ₀ = MAP.GenPath(x₀,xₑ,N);
plt = plot(x=φ₀[:,1],y=φ₀[:,2], Geom.point)
draw(SVG("InitialPath.svg"),plt);

grad = zeros(φ₀);
MAP.dS!(grad,φ₀,f,g,x₀,xₑ,N,n,Tspan/N)
plt = plot(x=grad[:,1],y=grad[:,2], Geom.point)
draw(SVG("InitialGrad.svg"),plt);

# Use Optim to optimise over path φ
resObj = MAP.MAP_Opt(f, g, x₀, xₑ, Tspan,N);
res = Optim.minimizer(resObj)
plt = plot(x=res[:,1],y=res[:,2], Geom.point)
draw(SVG("FinalPath.svg"),plt);

val = Optim.minimum(resObj)

Tvec = [];    Svec = [];
MAP.T_Opt!(Tvec,Svec, f,g,x₀,xₑ,N)
plt = plot(x=Tvec,y=Svec, Geom.point, Scale.y_log)
draw(SVG("TimeOptimisation.svg"),plt);
