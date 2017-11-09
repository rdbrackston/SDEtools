using Optim
# using ForwardDiff
using ReverseDiff
using Gadfly

# Specify the SDE

α = 0.5;
λ = 0.5α;
β = 0.0;
σ = 0.2;

f(x) = [2α*x[1]-4λ*x[1]^3-β; -4λ*x[2]^3];
g(x) = σ*[1.0; 1.0];

N = 101;    # Number of discrete points
T = 100.0;  # Time span
dT = T/N;
x₀ = [-1.0; 0.0];  # Start point
xₑ = [1.0; 0.0];  # End point
# xᵢ = [0.0; 0.0];  # Intermediate point

# Calculate the discretised piecewise path between the points
φ₀ = vcat(collect(linspace(x₀[1],xₑ[1],N))',collect(linspace(x₀[2],xₑ[2],N))');
φ₀ = φ₀[:,2:end-1];

# Define the action functional
function S(φ::AbstractArray)

    V = hcat(x₀,φ,xₑ);
    # print(size(V))
    act = 0.0;

    for iT = 1:N-1

        xT = 0.5*(V[:,iT+1]-V[:,iT]);
        fT = f(xT);
        gT = g(xT);

        for ii = 1:size(V)[1]
            dV = (V[ii,iT+1]-V[ii,iT])/dT;
            act += dT*(dV-fT[ii])^2/gT[ii]^2;
        end
    end

    return 0.5*act

end

store = zeros(φ₀);
# dS!(store,φ) = ForwardDiff.gradient!(store, S, φ);
dS!(store,φ) = ReverseDiff.gradient!(store, S, φ);

# Use Optim to optimise over path φ
resObj = Optim.optimize(S, dS!, φ₀, LBFGS());
res = Optim.minimizer(resObj);
