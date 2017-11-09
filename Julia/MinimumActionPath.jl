using Optim
using ForwardDiff
# using ReverseDiff
using Gadfly

# Specify the SDE

α = 0.5;
λ = 0.5α;
β = 0.0;
σ = 0.2;

f(x) = [2α*x[1]-4λ*x[1]^3-β; -4λ*x[2]^3];
g(x) = σ*[1.0; 1.0];

T = 100.0;  # Time span
x₀ = [-1.2; 0.0];  # Start point
xₑ = [1.2; 0.0];  # End point
# xᵢ = [0.0; 0.0];  # Intermediate point

function MAP(f::Function, g::Function, x₀::Vector, xₑ::Vector, T::Real)

    N = 101;    # Number of discrete points
    dT = T/N;

    # Calculate the discretised piecewise path between the points
    φ₀ = vcat(collect(linspace(x₀[1],xₑ[1],N))',collect(linspace(x₀[2],xₑ[2],N))');
    φ₀ = φ₀[:,2:end-1];

    # Define the action functional
    function S(φ::AbstractArray)

        V = hcat(x₀,φ,xₑ);
        # V = φ;

        action = 0.0;

        for iT = 1:N-1

            xT = 0.5*(V[:,iT+1]-V[:,iT]);
            fT = f(xT);
            gT = g(xT);

            for ii = 1:size(V)[1]
                dV = (V[ii,iT+1]-V[ii,iT])/dT;
                action += dT*(dV-fT[ii])^2/gT[ii]^2;
            end
        end

        return 0.5*action

    end

    store = zeros(φ₀);
    dS!(store,φ) = ForwardDiff.gradient!(store, S, φ);

    return Optim.optimize(S, dS!, φ₀, LBFGS());

end



# dS!(store,φ) = ReverseDiff.gradient!(store, S, φ);

# Use Optim to optimise over path φ
resObj = MAP(f, g, x₀, xₑ, T);
res = Optim.minimizer(resObj)
val = Optim.minimum(resObj)
