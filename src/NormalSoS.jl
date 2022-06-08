module NormalSoS

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials, CSDP, LinearAlgebra, Statistics
using Plots
gr()
export normdecomp


"""
Performs a two stage optimisation to try and obtain a Lyapunov function
satisfying the normal decomposition.
"""
function normdecomp(f, x, SDPsolver=CSDP.Optimizer, nIters=1, basis=:extended,
                    o=2, filter=true,
                    V::Union{DynamicPolynomials.Polynomial,Symbol}=:auto)

    n = length(f)

    # The Lyapunov function V(x):
    if basis == :minimal
        Z = minimalbasis(f,x)
    elseif basis == :extended
        Z = extendedbasis(f,x)
    elseif basis == :monomial
        Z = monomials(x,0:o)
    else
        println("Invalid basis selection. Exiting.")
        return
    end
    print("Chosen basis as:", "\n")
    print(Z, "\n")

    # Apply matrix constraint, ∇U⋅g ≤ 0.
    I = NormalSoS.eye(x)

    # Perform the first optimisation
    if V==:auto
        # V = normopt1(f,x,Z,SDPsolver,o)
        V = normopt2(f,x,Z,SDPsolver)
    end
    if any(isnan.(MultivariatePolynomials.coefficients(V)))
        println("Optimisation failed, exiting.")
        return zeros(size(Z))'*Z
    end

    # Now iterate to improve
    for ii=1:nIters

        U = V

        m = SOSModel(SDPsolver)
        @variable m ϵ
        @variable m α
        @variable m V Poly(Z)

        # Positive definiteness constraint
        @constraint m V ≥ ϵ*sum(x.^2)

        # Apply matrix constraint, ∇U⋅g ≤ 0.
        Mv = [-dot(differentiate(V, x),f); differentiate(V,x)]
        for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
        @constraint(m,Mv in PSDCone()) # Mv positive definite

        # Wynn inequality constraint
        @constraint m dot(differentiate(V,x),f+2*differentiate(U,x)) ≥
            α*dot(differentiate(U,x),f) + (1+α)*sum(differentiate(U,x).^2)

        @constraint m α ≥ 0
        @constraint m ϵ ≥ 0
        @objective m Min α

        status = optimize!(m)
        V = value(V)

        # If V now contains any NaN, return U instead
        if any(isnan.(MultivariatePolynomials.coefficients(V)))
            return filterterms(U,filter)
        end

        if checknorm(f,V,x) < 1e-10
            return filterterms(V,filter)
        end

    end

    return filterterms(V,filter)

end


"""
Obtains a Lyapunov function coming close to orthogonality. Uses a simple
lower bounding polynomial and maximises a single MultivariatePolynomials.coefficient.
"""
function normopt1(f, x, basis, SDPsolver=CSDP.Optimizer, o=2)
    # Single ϵ used for lower bound

    n = length(f)

    m = SOSModel(SDPsolver)
    @variable m ϵ

    @variable m V Poly(basis)

    # Positive definiteness constraint
    @constraint m V ≥ ϵ*sum(x.^o)

    # Apply matrix constraint, ∇U⋅g ≤ 0.
    I = NormalSoS.eye(x)
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
    @constraint(m,Mv in PSDCone()) # Mv positive definite

    @objective m Max ϵ
    status = optimize!(m)

    @show status
    @show sum(x.^o)
    @show value(ϵ)
    return value(V)

end


"""
Obtains a Lyapunov function coming close to orthogonality. Computes a suitable
lower bounding polynomial and maximises the sum of the MultivariatePolynomials.coefficients.
"""
function normopt2(f, x, basis, SDPsolver=CSDP.Optimizer, nonneg=false)
    # Vector ϵ used for lower bound

    n = length(f)

    m = SOSModel(SDPsolver)

    @variable m V Poly(basis)

    # Specify the lower bounding polynomial, bnd
    o = zeros(Int,1,n)
    # o[1] = maximum([degree(Vi,x[1]) for Vi in basis])
    if n>1
        bnd = x[1]    # Initialise the type of bnd
        # First add the maximum even order polynomial for each x[ii]
        for ii=1:n
            o[ii] = maximum(filter(iseven,[degree(Vi,x[ii]) for Vi in basis]))
            if o[ii]>0;    bnd += x[ii]^o[ii];    end
        end
        # Now add any fully even terms of mixed x[ii], e.g. x[1]^2*x[2]^2
        for bi in basis
            e = exponents(bi)
            if all(y->y==e[1]&&e[1]%2==0, e) && e[1]!=0
                bnd += bi
            end
        end
        bnd -= x[1]
        b = length(bnd)
    else
        b = 1
        o[1] = maximum([degree(Vi,x[1]) for Vi in basis])
        bnd = x[1]^o[1]
        bnd = [bnd]
    end

    # Positive definitness constraint on V
    @variable m ϵ[1:b]
    for ii=1:b; @constraint m ϵ[ii] ≥ 0; end

    # Make Mᵥ for matrix constraint, ∇U⋅fᵥ ≤ 0.
    I = NormalSoS.eye(x)
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)]
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end

    # if nonneg
        # Generate the set for non-negative x
        # s = @set 0 ≤ x[1]
        # for ii=2:n
        #     s = @set 0 ≤ x[ii] && s;
        # end
        # @constraint(m, V≥sum([ϵ[ii]*bnd[ii] for ii=1:b]), domain=s)
        # @constraint(m, Mv in PSDCone(), domain=s); # Mv positive definite
    # else
        @constraint m V ≥ sum([ϵ[ii]*bnd[ii] for ii=1:b])
        @constraint(m,Mv in PSDCone()) # Mv positive definite
    # end

    @objective m Max sum(ϵ)
    # TT = STDOUT; # save original STDOUT stream
    # redirect_stdout()
    status = optimize!(m)
    # redirect_stdout(TT)

    @show status
    @show bnd
    @show value.(ϵ)
    return value.(V)

end


"""
A basic function to obtain a Lyapunov function for ODE f. The function can
optionally be applied only for positive x, e.g. for a system describing a
chemical reaction network.
"""
function lyapunov(f, x, SDPsolver=CSDP.Optimizer, o=2, nonneg=false)

    n = length(f)

    m = SOSModel(SDPsolver)
    # @variable m ϵ
    @variable m V Poly(monomials(x,0:o))

    # Positive definiteness constraint
    @constraint m V ≥ sum(x.^2)

    # Standard Lyapunov constraint
    P = dot(differentiate(V, x),f)

    # if nonneg
        # s = @set 0 ≤ x[1]
        # for ii=2:n
        #     s = @set 0 ≤ x[ii] && s
        # end
        # @constraint(m, P ≤ 0, domain=s)
    # else
        @constraint m P ≤ 0
    # end

    status = optimize!(m)
    @show status
    # @show value(ϵ)
    return value(V)

end


"""
A minimal working example to obtain a Lyapunov function. Originally made for the
question in julia discourse.
"""
function minlyapunov(f,x,o=2)

    n = length(x)

    m = SOSModel(CSDPSolver())
    @variable m ϵ
    @variable m V Poly(monomials(x,0:o))

    # Make the semialgebraicset of non-negative x
    # s = @set 0 ≤ x[1]
    # for ii=2:n
    #     s = @set 0 ≤ x[ii] && s
    # end

    # Positive definiteness constraint
    # @constraint(m, V ≥ ϵ*sum(x.^2), domain=s)
    @constraint(m, V ≥ ϵ*sum(x.^2))
    @constraint m ϵ ≥ 0

    # Apply matrix constraint, ∇U⋅fᵥ ≤ 0.
    I = NormalSoS.eye(x)
    Mv = [-dot(differentiate(V,x),f); differentiate(V,x)]
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end

    @constraint(m, Mv in PSDCone(), domain=s) # Mv positive definite

    status = optimize!(m)
    @show(status)

    return value(V)

end


"""
Obtain a minimal basis for V, based on the terms required to form f.
"""
function minimalbasis(f,x)

    basis = 1.0

    # Loop over the elements of f
    for (i,fi) in enumerate(f)
        fTmp = 0
        for elem in fi
            fTmp += MultivariatePolynomials.coefficient(elem)*elem
        end
        basis = basis + fTmp*x[i]
    end
    basis = monomials(basis)

    return basis

end


"""
Obtain an extended basis for V, starting from the terms required to form f, but
with the addition of cross-terms up to the maximum order in each xᵢ.
"""
function extendedbasis(f,x)

    # Start from the minimal basis
    basis = minimalbasis(f,x)
    n = length(x)

    # Find maximum total degree
    d = maximum([degree(Vi) for Vi in basis])

    # Now find maximum individual degree for each xᵢ, o[ii]
    o = zeros(Int,1,n)
    for ii=1:n
        o[ii] = maximum([degree(Vi,x[ii]) for Vi in basis])
    end

    # Add all mixed terms up to total degree d, and with maximum individual degree o[ii]
    basis = monomials(x,0:d, m->check(m,x,o))
    return basis

end


"""
Define the function that confirms that all for monomial term m, each variable
x[ii] has degree less than o[ii].
"""
function check(m,x,o)
    for ii=1:length(x)
        if degree(m,x[ii])>o[ii]
            return false
        end
    end
    return true
end


"""
Return the identity matrix of the same dimension and Type as x.
"""
function eye(x)
    n = length(x)
    v = [1.0+0.0x[1]]
    for i = 1:n-1
        push!(v,1.0+0.0x[1])
    end
    return diagm(v)
end

"""
Function to implement the new iterative upper bound method
"""
function upperbound(f, x, linear=true, niters=2)

    n = length(f)
    q = 2

    # The chosen basis for V(x):
    Z = extendedbasis(f,x)
    print("Chosen basis as:", "\n")
    print(Z, "\n")

    # Correctly typed identity matrix
    I = eye(x)

    # Function c required in optimization
    c(g,h) = dot(differentiate(h,x), -f-2*differentiate(g,x)) + sum(differentiate(g,x).^2)

    # Obtain initial guess for V₀
    V = initialguess(f,x,linear)

    # Begin iterative loop
    for ii=1:niters

        m = SOSModel(CSDP.Optimizer)
        @variable m W Poly(Z)

        @variable m ϵ[1:n]
        for ii=1:n; @constraint m -ϵ[ii] ≥ 0; end

        @constraint m c(V,W) - sum([ϵ[ii]*x[ii]^2 for ii=1:n]) ≥ 0
        @constraint m c(V,W) - c(V,V) ≥ 0
        @constraint m W ≥ 0
        @constraint m -c(V,W) ≥ 0

        @objective m Max sum(ϵ)
        status = optimize!(m)
        @show(status)

        # Set new V₀ as the optimization output
        V = value(W)

    end

    return V

end


"""
Function to provide the initial guess V₀ for the upper bound method
"""
function initialguess(f, x, linear=true)

    if linear
        println("Obtaining initial guess under condition of linearity.")

        m = SOSModel(CSDP.Optimizer)
        # @variable m α
        @variable m β
        α = 1.0e0

        Exp = -β*α*sum(x.^2) + 0.25*sum((α*x-f).^2)
        # @show α*x-f
        @constraint m Exp ≤ 0
        @constraint m β ≥ 0
        # @constraint m α ≥ 0

        @objective m Min β
        status = optimize!(m)
        @show(status)

        # @show value(α)
        return value(β)*sum(x.^2)

    else
        println("Initial program for nonlinear case.")

        d = maximum([maximum([degree(t) for t in fi]) for fi in f])
        n = length(f)
        u = 0.5*(x*sum(x.^(d-1)) - f)
        @show(u)
        Z = extendedbasis(f,x)

        m = SOSModel(CSDP.Optimizer)
        @variable m V Poly(Z)
        @variable m ϵ[1:n]

        I = NormalSoS.eye(x)
        M = [-dot(differentiate(V,x),f+2*u) - sum([ϵ[ii]*x[ii]^2 for ii=1:n]);
              u]
        for ii=1:n; M = hcat(Mv, [u[ii]; I[:,ii]]); end
        @constraint(m,M in PSDCone())

        @objective m Min sum(ϵ)
        status = optimize!(m)
        @show(status)

        return
    end

end


"""
Plot the first two dimensions of the landscape U. A quiver plot of the vector
field f may optionally be included.
"""
function plotlandscape(f, U, x, lims, vectors=false, scl=0.05)

    Ng = 100
    Nds = 10

    # Assemble the grid of points
    xv = collect(LinRange(lims[1][1],lims[1][2],Ng))
    yv = collect(LinRange(lims[2][1],lims[2][2],Ng))

    # Evaluate U using an array comprehension then plot
    for ii=1:length(x)-2
        U = subs(U,x[ii+2]=>0.0)
        f = subs(f,x[ii+2]=>0.0)
    end
    Umat = [coefficients(subs(U, x[1]=>xv[ii], x[2]=>yv[jj]))[1] for ii=1:Ng, jj=1:Ng]
    plt = Plots.contour(xv,yv,Umat'.-minimum(Umat), aspect_ratio=:equal,
                xlims=(lims[1][1],lims[1][2]), ylims=(lims[2][1],lims[2][2]))

    # If desired, evaluate f using an array comprehension then plot
    gU = differentiate(U,x) # Gradient
    fU = f + gU             # Curl component
    if vectors
        xm = vec([xv[ii] for jj=1:Nds:Ng, ii=1:Nds:Ng])
        ym = vec([yv[ii] for ii=1:Nds:Ng, jj=1:Nds:Ng])
        fMat = vec([(scl.*coefficients(subs(f[1], x[1]=>xv[ii], x[2]=>yv[jj]))[1],
                     scl.*coefficients(subs(f[2], x[1]=>xv[ii], x[2]=>yv[jj]))[1])
               for jj=1:Nds:Ng, ii=1:Nds:Ng])
        gMat = vec([(scl.*coefficients(subs(-gU[1], x[1]=>xv[ii], x[2]=>yv[jj]))[1],
                     scl.*coefficients(subs(-gU[2], x[1]=>xv[ii], x[2]=>yv[jj]))[1])
               for jj=1:Nds:Ng, ii=1:Nds:Ng])
        cMat = vec([(scl.*coefficients(subs(fU[1], x[1]=>xv[ii], x[2]=>yv[jj]))[1],
                     scl.*coefficients(subs(fU[2], x[1]=>xv[ii], x[2]=>yv[jj]))[1])
               for jj=1:Nds:Ng, ii=1:Nds:Ng])
        Plots.quiver!(xm,ym, quiver=gMat, color=:black)
        Plots.quiver!(xm,ym, quiver=cMat, color=:black)
    end

    return plt

end


"""
Generate a quiver plot of the vector field f.
"""
function plotvectors(f, x, lims, scl=0.05)

    Ng = 100
    Nds = 10

    # Assemble the grid of points
    xv = collect(LinRange(lims[1][1],lims[1][2],Ng))
    yv = collect(LinRange(lims[2][1],lims[2][2],Ng))

    xm = vec([xv[ii] for jj=1:Nds:Ng, ii=1:Nds:Ng])
    ym = vec([yv[ii] for ii=1:Nds:Ng, jj=1:Nds:Ng])
    fMat = vec([(scl.*Float64(subs(f[1], x[1]=>xv[ii], x[2]=>yv[jj])),
                 scl.*Float64(subs(f[2], x[1]=>xv[ii], x[2]=>yv[jj])))
           for jj=1:Nds:Ng, ii=1:Nds:Ng])
    plt = Plots.quiver(xm,ym, quiver=fMat, color=:black, aspect_ratio=:equal,
                xlims=(lims[1][1],lims[1][2]), ylims=(lims[2][1],lims[2][2]))

    return plt

end


"""
Generate a metric that quantifies the orthogonality of fᵤ with ∇U.
"""
function checknorm(f, U, x)

    # Evaluate ∇U
    ∇U = differentiate(U,x)
    res = dot(∇U,f) + dot(∇U,∇U)

    return mean(abs.(MultivariatePolynomials.coefficients(res))) /
            mean([mean(abs.(MultivariatePolynomials.coefficients(f[ii]))) for ii=1:length(x)])
end


"""
Function to remove terms with very small MultivariatePolynomials.coefficients
"""
function filterterms(U, filter=true, tol=1e-4)

    if filter
        U2 = 1

        idxs = abs.(MultivariatePolynomials.coefficients(U)).>tol
        for (ii,trm) in enumerate(U)
            if idxs[ii]
                U2 += trm
            end
        end

        U2 -= 1
        return U2
    else
        return U
    end

end


end
