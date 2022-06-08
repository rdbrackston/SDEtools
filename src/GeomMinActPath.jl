module gMAM
# Module for computing the geometric minimum action path.
# Currently only computation of the action is working.

using Optim, Printf
using ForwardDiff


"""
Evaluate the time-independent geometric action for a discretised path φ.
"""
function action(φ::AbstractArray,
           f::Function, g::Function)

    # Loop over the elements in φ
    V = φ
    (N,n) = size(V)

    action = 0.0
    for jj = 1:N-1

        dl² = 0.0
        b = 0.0
        c = 0.0
        d = 0.0

        for ii = 1:n
            dl² += (V[jj+1,ii]-V[jj,ii])^2

            b += ((f(V[jj+1,:])[ii]+f(V[jj,:])[ii])/
                                        (g(V[jj+1,:])[ii]+g(V[jj,:])[ii]))^2

            c += (V[jj+1,ii]-V[jj,ii])*(f(V[jj+1,:])[ii]+f(V[jj,:])[ii])/
                                        (g(V[jj+1,:])[ii]+g(V[jj,:])[ii])^2

            d += (V[jj+1,ii]-V[jj,ii])^2/(g(V[jj+1,:])[ii]+g(V[jj,:])[ii])^2
        end

        action += (sqrt(b)-c/sqrt(d))*sqrt(dl²)
    end

    return 0.5*action

end


"""
Compute the gradient of the geometric action using automatic differentiation.
"""
function actiongradient!(store::AbstractArray, φ::AbstractArray,
           f::Function, g::Function,
           X₀::AbstractArray, Xₑ::AbstractArray,
           N::Int64, n::Int64)

    func = v-> action(v,f,g,X₀,Xₑ,N,n)

    return ForwardDiff.gradient(func,φ)

end


"""
Perform the gradient-based optimisation over the discretised path.
"""
function optimalpath(f::Function, g::Function,
                     x₀::AbstractArray, xₑ::AbstractArray, N::Signed,
                     φ₀::Union{AbstractArray,Symbol}=:auto, reRuns::Signed=1)

    n = length(x₀)       # Number of state dimensions
    # method = LBFGS()
    method = BFGS()

    # Evaluate the initial path if not given
    if φ₀==:auto
        φ₀ = makepath(x₀,xₑ,N)
    end

    # Define the anonymous action function
    S_opt = φ->action(φ, f,g,x₀,xₑ,N,n)

    store = zeros(φ₀)
    dS_opt! = (store,φ)->actiongradient!(store,φ, f,g,x₀,xₑ,N,n)

    # Perform the optimisation and print some info
    OptStruct = Optim.optimize(S_opt, dS_opt!, φ₀, method)
    println(@printf("Optimisation gives S=%.2f",Optim.minimum(OptStruct)))

    ii = 0
    while !Optim.converged(OptStruct) && ii<reRuns
        println("Optimisation is not converged, rerunning...")
        φ₀ = Optim.minimizer(OptStruct)
        OptStruct = Optim.optimize(S_opt, dS_opt!, φ₀, method)
        println(@printf("Optimisation gives S=%.2f",Optim.minimum(OptStruct)))
        ii += 1
    end

    if Optim.converged(OptStruct)
        println("Optimisation is converged.")
    else
        println("Optimisation not converged")
    end
    return OptStruct
end


end # module
