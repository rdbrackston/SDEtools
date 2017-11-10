using Optim
using ForwardDiff
# using ReverseDiff
using Gadfly

# Specify the SDE

α = 0.5;
λ = 0.5α;
β = 0.0;
σ = 0.2;

# f(x) = [2α*x[1]-4λ*x[1]^3-β; -4λ*x[2]^3];
f(x::Vector) = [2α*x[1]-4λ*x[1]^3-β];
# g(x) = σ*[1.0; 1.0];
g(x::Vector) = [σ];

Tspan = 100.0;  # Time span
x₀ = [-1.0];  # Start point (row vector)
xₑ = [1.0];  # End point (row vector)
# xᵢ = [0.0, 0.0];  # Intermediate point

# Function to generate the initial path
function GenPath(X₀::AbstractArray, Xₑ::AbstractArray, N::Int64)

    n = length(X₀);

    # Calculate the discretised piecewise path between the points
    φ₀ = zeros(Float64, N,n);
    for ii=1:n
        φ₀[:,ii] = collect(linspace(x₀[ii],xₑ[ii],N));
    end
    if n==1
        φ₀ = vec(φ₀[2:end-1,:]);    # Subset of modifiable points as column vector
    else
        φ₀ = φ₀[2:end-1,:];    # Subset of modifiable points as column vector
    end

    return φ₀
end

# Define the action functional
function S(φ::AbstractArray,
           X₀::AbstractArray,
           Xₑ::AbstractArray,
           N::Int64,
           n::Int64,
           dt::Float64)

    # V = hcat(x₀,φ,xₑ);
    V = vcat(X₀,φ,Xₑ);

    action = 0.0;
    for iT = 1:N-1

        xT = 0.5*(V[iT+1,:]-V[iT,:]);
        fT = f(xT);
        DT = g(xT).^2;

        for ii = 1:n
            dV = (V[iT+1,ii]-V[iT,ii])/dt;
            action += dt*(dV-fT[ii])^2/DT[ii];
        end
    end

    return 0.5*action
end

# Analytical gradient of action functional
function dS!(store::AbstractArray,
             φ::AbstractArray,
             X₀::AbstractArray,
             Xₑ::AbstractArray,
             N::Int64,
             n::Int64,
             dt::Float64)

    V = vcat(x₀,φ,xₑ);    # Complete vector of points

    # Evaluate for first segment
    xT = 0.5*(V[2,:]-V[1,:]);
    fT_prev = f(xT);
    DT_prev = g(xT).^2;
    dφ_prev = (V[2,:]-V[1,:])/dt;

    # Form the arrays of functions
    F = Array{Function}(n);
    D = Array{Function}(n);
    df_prev = zeros(Float64, n,n);
    dD_prev = zeros(Float64, n,n);
    df_next = zeros(Float64, n,n);
    dD_next = zeros(Float64, n,n);
    for ii = 1:n
        F[ii] = x->f(x)[ii];
        D[ii] = x->(g(x)[ii])^2;

        df_prev[:,ii] = ForwardDiff.gradient(F[ii],xT);
        dD_prev[:,ii] = ForwardDiff.gradient(D[ii],xT);
    end

    # Loop over the items in store
    for ik = 1:N-2

        # Evaluate the quantities after this point
        xT = 0.5*(V[ik+2,:]-V[ik+1,:]);
        fT_next = f(xT);
        DT_next = g(xT).^2;
        dφ_next = (V[ik+2,:]-V[ik+1,:])/dt;
        for ii = 1:n
            df_next[:,ii] = ForwardDiff.gradient(F[ii],xT);
            dD_next[:,ii] = ForwardDiff.gradient(D[ii],xT);
        end

        for ii = 1:n    # Loop over coordinate dimensions
            store[ik,ii] = 0.0;

            store[ik,ii] += (dφ_prev[ii] - fT_prev[ii])/DT_prev[ii];
            store[ik,ii] -= (dφ_next[ii] - fT_next[ii])/DT_next[ii];

            for ij = 1:n
                store[ik,ii] -= 0.5dt*(dφ_next[ij]-fT_next[ij])*df_next[ii,ij]/DT_next[ij];
                store[ik,ii] -= 0.5dt*(dφ_prev[ij]-fT_prev[ij])*df_prev[ii,ij]/DT_prev[ij];
                store[ik,ii] -= 0.25dt*dD_next[ii,ij]*(dφ_next[ij]-fT_next[ij])^2/DT_next[ij]^2;
                store[ik,ii] -= 0.25dt*dD_prev[ii,ij]*(dφ_prev[ij]-fT_prev[ij])^2/DT_prev[ij]^2;
            end
        end

        fT_prev = deepcopy(fT_next);
        DT_prev = deepcopy(DT_next);
        dφ_prev = deepcopy(dφ_next);
        df_prev = deepcopy(df_next);
        dD_prev = deepcopy(dD_next);
        df_prev = deepcopy(df_next);
        dD_prev = deepcopy(dD_next);
    end
end

function MAP(f::Function, g::Function, x₀::Vector, xₑ::Vector, τ::Real)

    N = 101;    # Number of discrete points including start and end
    dτ = τ/N;
    n = length(x₀);       # Number of state dimensions

    φ₀ = GenPath(x₀,xₑ,N);

    # Define the anonymous action function
    S_opt = φ->S(φ, x₀,xₑ,N,n,dτ);

    store = zeros(φ₀);
    # dS!(store,φ) = ForwardDiff.gradient!(store, S, φ);
    dS_opt! = (store, φ)->dS!(store,φ, x₀,xₑ,N,n,dτ);

    return Optim.optimize(S_opt, dS_opt!, φ₀, LBFGS());
    # return Optim.optimize(S, φ₀, LBFGS());
end

# dS!(store,φ) = ReverseDiff.gradient!(store, S, φ);

# Use Optim to optimise over path φ
resObj = MAP(f, g, x₀, xₑ, Tspan);
res = Optim.minimizer(resObj)
val = Optim.minimum(resObj)
