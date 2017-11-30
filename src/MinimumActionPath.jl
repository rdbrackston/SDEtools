module MAP

using Optim
using ForwardDiff

abstract type Tracker end

struct OptimisationTracker <: Tracker
    pointVec::AbstractArray{AbstractFloat};
    valueVec::AbstractArray{AbstractFloat};
    paths::AbstractArray{AbstractArray};
    convergences::AbstractArray{Bool};
end

function updatetracker!(trackerObj::OptimisationTracker,
                        tVal::AbstractFloat, optStruct)

    push!(trackerObj.pointVec,tVal);
    push!(trackerObj.valueVec,Optim.minimum(optStruct));
    push!(trackerObj.paths,Optim.minimizer(optStruct));
    push!(trackerObj.convergences,Optim.converged(optStruct));
end

# Function to generate the initial path
function makepath(X₀::AbstractArray, Xₑ::AbstractArray, N::Int64)

    n = length(X₀);

    # Calculate the discretised piecewise path between the points
    φ₀ = zeros(Float64, N,n);
    for ii=1:n
        φ₀[:,ii] = collect(linspace(X₀[ii],Xₑ[ii],N));
    end
    if n==1
        φ₀ = vec(φ₀[2:end-1,:]);    # Subset of modifiable points as column vector
    else
        φ₀ = φ₀[2:end-1,:];    # Subset of modifiable points as column vector
    end

    return φ₀
end

# Define the action functional
function action(φ::AbstractArray,
           f::Function, g::Function,
           X₀::AbstractArray, Xₑ::AbstractArray,
           N::Int64, n::Int64, dt::Float64)

    # V = hcat(x₀,φ,xₑ);
    V = vcat(X₀,φ,Xₑ);

    action = 0.0;
    for iT = 1:N-1

        xT = 0.5*(V[iT+1,:]+V[iT,:]);
        fT = f(xT);
        DT = g(xT).^2;

        for ii = 1:n
            dV = (V[iT+1,ii]-V[iT,ii])/dt;
            action += dt*(dV-fT[ii])^2/DT[ii];
        end
    end

    return 0.5*action
end

# Analytical gradient of the action functional
function actiongradient!(store::AbstractArray, φ::AbstractArray,
                         f::Function, g::Function,
                         X₀::AbstractArray, Xₑ::AbstractArray,
                         N::Int64, n::Int64, dt::Float64)

    V = vcat(X₀,φ,Xₑ);    # Complete vector of points

    # Evaluate for first segment
    xT = 0.5*(V[2,:]+V[1,:]);
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
        xT = 0.5*(V[ik+2,:]+V[ik+1,:]);
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

function optimalpath(f::Function, g::Function,
                     x₀::AbstractArray, xₑ::AbstractArray,
                     τ::Real, N::Signed,
                     φ₀::Union{AbstractArray,Symbol}=:auto, reRuns::Signed=1)

    dτ = τ/N;
    n = length(x₀);       # Number of state dimensions
    # method = ConjugateGradOptim.convergedient();
    method = LBFGS();
    # method = BFGS();

    # Evaluate the initial path if not given
    if φ₀==:auto
        φ₀ = makepath(x₀,xₑ,N);
    end

    # Define the anonymous action function
    S_opt = φ->action(φ, f,g,x₀,xₑ,N,n,dτ);

    store = zeros(φ₀);
    dS_opt! = (store,φ)->actiongradient!(store,φ, f,g,x₀,xₑ,N,n,dτ);

    # Perform the optimisation and print some info
    OptStruct = Optim.optimize(S_opt, dS_opt!, φ₀, method);
    println(@sprintf("Optimisation for T=%.2f gives S=%.2f",τ,Optim.minimum(OptStruct)));

    ii = 0;
    while !Optim.converged(OptStruct) && ii<reRuns
        println("Optimisation is not converged, rerunning...")
        φ₀ = Optim.minimizer(OptStruct)
        OptStruct = Optim.optimize(S_opt, dS_opt!, φ₀, method);
        println(@sprintf("Optimisation for T=%.2f gives S=%.2f",τ,Optim.minimum(OptStruct)));
        ii += 1;
    end

    if Optim.converged(OptStruct)
        println("Optimisation is converged.")
    else
        println("Optimisation not converged")
    end
    return OptStruct
end

function optimaltime!(trackerObj::OptimisationTracker,
                      f::Function, g::Function,
                      x₀::AbstractArray, xₑ::AbstractArray,
                      TBounds::NTuple{2,Real}, nPoints::Signed)
    # Function to optimise over the time duration using a golden section line
    # search algorithm

    nIter = 30;
    ϕ = 0.5(1+sqrt(5)); # Golden ratio
    ii = 0;

    # Define the recursive golden section algorithm
    function goldsectsearch(a::Float64,b::Float64,c::Float64,
                            Sa::Float64,Sb::Float64,Sc::Float64,
                            φa::AbstractArray, φb::AbstractArray,
                            tol::Float64)

        if (c-b)>(b-a)
            # d = b + (2-ϕ)*(c-b);
            d = (c+ϕ*b)/(1+ϕ);
        else
            d = (b+ϕ*a)/(1+ϕ);
        end

        if ii==nIter
            # Return the optimisasion result half way through the interval.
            println("Evaluating final point")
            OptStruct = optimalpath(f,g,x₀,xₑ,0.5(a+c),nPoints,φa);

            updatetracker!(trackerObj, 0.5(a+c),OptStruct)
            # push!(trackerObj.pointVec,0.5(a+c));
            # push!(trackerObj.valueVec,Optim.minimum(OptStruct));
            return OptStruct
        end

        ii += 1;
        println(@sprintf("Starting iteration no. %i out of %i", ii,nIter));
        # Evaluate the optimisation at the new point, using previous result at
        # closest existing point as the initial guess
        if (c-b)>(b-a)
            tmp = optimalpath(f,g,x₀,xₑ,d,nPoints,φb);
        else
            tmp = optimalpath(f,g,x₀,xₑ,d,nPoints,φa);
        end
        φd = Optim.minimizer(tmp);
        Sd = Optim.minimum(tmp);

        updatetracker!(trackerObj, d,tmp)
        # push!(trackerObj.pointVec,d);
        # push!(trackerObj.valueVec,Sd);

        if d>b
            if Sd>Sb
                return goldsectsearch(a,b,d, Sa,Sb,Sd, φa,φb, tol)
            else
                return goldsectsearch(b,d,c, Sb,Sd,Sc, φb,φd, tol)
            end
        else
            if Sd>Sb
                return goldsectsearch(d,b,c, Sd,Sb,Sc, φd,φb, tol)
            else
                return goldsectsearch(a,d,b, Sa,Sd,Sb, φa,φd, tol)
            end
        end
    end

    println("Evaluating initial triplet...")
    τL = TBounds[1];    τU = TBounds[2]; # Initial bounds
    τM = (τU+ϕ*τL)/(1+ϕ);    # Third point in starting triple

    tmp = optimalpath(f,g,x₀,xₑ,τL,nPoints);
    fL = Optim.minimum(tmp);
    φL = Optim.minimizer(tmp);
    updatetracker!(trackerObj, τL,tmp);

    tmp = optimalpath(f,g,x₀,xₑ,τM,nPoints);
    fM = Optim.minimum(tmp);
    φM = Optim.minimizer(tmp);
    updatetracker!(trackerObj, τM,tmp);

    tmp = optimalpath(f,g,x₀,xₑ,τU,nPoints);
    fU = Optim.minimum(tmp);
    updatetracker!(trackerObj, τU,tmp);

    # push!(trackerObj.pointVec,τL); push!(trackerObj.pointVec,τM); push!(trackerObj.pointVec,τU);
    # push!(trackerObj.valueVec,fL); push!(trackerObj.valueVec,fM); push!(trackerObj.valueVec,fU);
    return goldsectsearch(τL,τM,τU, fL,fM,fU, φL,φM, 0.0)

end

end
