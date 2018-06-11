using Optim
using Plots
gr()

include("../src/MinimumActionPath.jl");    using MAP

α = 0.5;
λ = 0.5α;
β = -0.05;
σ = 0.4;

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

Tspan = 31.0;  # Time span
N = 100;

# Examine the path
φ₀ = MAP.makepath(x₀,xₑ,N);
plot(φ₀[:,1],φ₀[:,2], line=0,marker=3)

grad = zeros(φ₀);
MAP.dS!(grad,φ₀,f,g,x₀,xₑ,N,n,Tspan/N)
plot(grad[:,1],grad[:,2])

# Use Optim to optimise over path φ
resObj = MAP.optimalpath(f, g, x₀, xₑ, Tspan,N);
res = Optim.minimizer(resObj)

Tvec = [];    Svec = [];
resObj = MAP.optimaltime!(Tvec,Svec, f,g,x₀,xₑ,(1.,100.),N)
res = Optim.minimizer(resObj)
plot(Tvec,Svec)

plot(res[:,1],res[:,2], line=0,marker=3)

resObjrev = MAP.optimaltime!(Tvec,Svec, f,g,xₑ,x₀,(5.,100.),N)
resrev = Optim.minimizer(resObj)
plot(resrev[:,1],resrev[:,2], line=0,marker=3)
