module NormalSoS

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials, CSDP
using Plots
gr()
export normdecomp

# To do:
#  - Add plotting functions
#  - Add function to check normality


"""
Performs a two stage optimisation to try and obtain a Lyapunov function
satisfying the normal decomposition.
"""
function normdecomp(f, x, SDPsolver=CSDPSolver(), nIters=1, o=2, basis=:minimal,
                    V::Union{DynamicPolynomials.Polynomial,Symbol}=:auto)

    n = length(f);

    # The Lyapunov function V(x):
    if basis == :minimal
        Z = minimalbasis(f,x);
    elseif basis == :monomial
        Z = monomials(x,0:o);
    else
        println("Invalid basis selection. Exiting.")
        return
    end
    print("Chosen basis as:", "\n")
    print(Z, "\n")

    # Apply matrix constraint, ∇U⋅g ≤ 0.
    I = NormalSoS.eye(x);

    # Perform the first optimisation
    if V==:auto
        # V = normopt1(f,x,Z,SDPsolver,o);
        V = normopt2(f,x,Z,SDPsolver);
    end

    # Now iterate to improve
    for ii=1:nIters

        U = V;

        m = SOSModel(solver=SDPsolver);
        @variable m ϵ
        @variable m α
        @polyvariable m V Z

        # Positive definiteness constraint
        @polyconstraint m V ≥ ϵ*sum(x.^2);

        # Apply matrix constraint, ∇U⋅g ≤ 0.
        Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
        for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
        @SDconstraint m Mv ⪰ 0 # Mv positive definite

        # Wynn inequality constraint
        @polyconstraint m dot(differentiate(V,x),f+2*differentiate(U,x)) ≥
            α*dot(differentiate(U,x),f) + (1+α)*sum(differentiate(U,x).^2)

        @constraint m α ≥ 0
        @constraint m ϵ ≥ 0
        @objective m Min α

        status = solve(m);
        V = getvalue(V);

    end

    return V

end


"""
Obtains a Lyapunov function coming close to orthogonality. Uses a simple
lower bounding polynomial and maximises a single coefficient.
"""
function normopt1(f, x, basis, SDPsolver=CSDPSolver(), o=2)
    # Single ϵ used for lower bound

    n = length(f);

    m = SOSModel(solver=SDPsolver);
    @variable m ϵ

    @polyvariable m V basis

    # Positive definiteness constraint
    @polyconstraint m V ≥ ϵ*sum(x.^o);

    # Apply matrix constraint, ∇U⋅g ≤ 0.
    I = NormalSoS.eye(x);
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
    @SDconstraint m Mv ⪰ 0 # Mv positive definite

    @objective m Max ϵ
    status = solve(m);

    @show status
    @show sum(x.^o)
    @show getvalue(ϵ)
    return getvalue(V)

end


"""
Obtains a Lyapunov function coming close to orthogonality. Computes a suitable
lower bounding polynomial and maximises the sum of the coefficients.
"""
function normopt2(f, x, basis, SDPsolver=CSDPSolver(), nonneg=false)
    # Vector ϵ used for lower bound

    n = length(f);

    m = SOSModel(solver=SDPsolver);

    @polyvariable m V basis

    # Specify the lower bounding polynomial, bnd
    o = zeros(Int,1,n)
    # o[1] = maximum([degree(Vi,x[1]) for Vi in basis])
    if n>1
        bnd = x[1];    # Initialise the type of bnd
        # First add the maximum even order polynomial for each x[ii]
        for ii=1:n;
            o[ii] = maximum(filter(iseven,[degree(Vi,x[ii]) for Vi in basis]));
            if o[ii]>0;    bnd += x[ii]^o[ii];    end
        end
        # Now add any fully even terms of mixed x[ii], e.g. x[1]^2*x[2]^2
        for bi in basis
            e = exponents(bi);
            if all(y->y==e[1]&&e[1]%2==0, e) && e[1]!=0
                bnd += bi
            end
        end
        bnd -= x[1];
        b = length(bnd);
    else
        b = 1;
        o[1] = maximum([degree(Vi,x[1]) for Vi in basis])
        bnd = x[1]^o[1];
        bnd = [bnd];
    end

    # Positive definitness constraint on V
    @variable m ϵ[1:b];
    for ii=1:b; @constraint m ϵ[ii] ≥ 0; end

    # Make Mᵥ for matrix constraint, ∇U⋅fᵥ ≤ 0.
    I = NormalSoS.eye(x);
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end

    if nonneg
        # Generate the set for non-negative x
        s = @set 0 ≤ x[1]
        for ii=2:n
            s = @set 0 ≤ x[ii] && s;
        end
        @polyconstraint(m, V≥sum([ϵ[ii]*bnd[ii] for ii=1:b]), domain=s)
        @constraint(m, Mv in PSDCone(), domain=s); # Mv positive definite
    else
        @polyconstraint m V ≥ sum([ϵ[ii]*bnd[ii] for ii=1:b])
        @SDconstraint m Mv ⪰ 0 # Mv positive definite
    end

    @objective m Max sum(ϵ)
    status = solve(m);

    @show status
    @show bnd
    @show getvalue(ϵ)
    return getvalue(V)

end


"""
A basic function to obtain a Lyapunov function for ODE f. The function can be
applied only for positive x, e.g. for a system describing a chemical reaction
network.
"""
function lyapunov(f, x, SDPsolver=CSDPSolver(), o=2, nonneg=false)

    n = length(f);

    m = SOSModel(solver=SDPsolver);
    # @variable m ϵ
    @polyvariable m V monomials(x,o);

    # Positive definiteness constraint
    @polyconstraint m V ≥ sum(x.^2);
    # @constraint m ϵ ≥ 0
    # @objective m Max ϵ

    # Standard Lyapunov constraint
    P = dot(differentiate(V, x),f);

    if nonneg
        s = @set 0 ≤ x[1]
        for ii=2:n
            s = @set 0 ≤ x[ii] && s;
        end
        @constraint(m, P ≤ 0, domain=s)
    else
        @constraint m P ≤ 0;
    end

    status = solve(m)
    @show status
    # @show getvalue(ϵ)
    return getvalue(V)

end


"""
A minimal working example to obtain a Lyapunov function. Originally made for the
question in julia discourse.
"""
function minlyapunov(f,x)

    n = length(x);

    m = SOSModel(solver=CSDPSolver());
    @variable m ϵ
    @polyvariable m V monomials(x,2);

    # Make the semialgebraicset of non-negative x
    s = @set 0 ≤ x[1]
    for ii=2:n
        s = @set 0 ≤ x[ii] && s;
    end

    # Positive definiteness constraint
    @polyconstraint(m, V ≥ ϵ*sum(x.^2), domain=s);
    @constraint m ϵ ≥ 0

    # Apply matrix constraint, ∇U⋅fᵥ ≤ 0.
    I = NormalSoS.eye(x);
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end

    @constraint(m, Mv in PSDCone(), domain=s); # Mv positive definite

    status = solve(m);
    @show(status)

    return getvalue(V)

end


"""
Obtain a minimal basis for V, based on the terms required to form f.
"""
function minimalbasis(f,x)

    basis = 1.0;

    # Loop over the elements of f
    for (i,fi) in enumerate(f)
        fTmp = [];
        for elem in fi
            fTmp += coefficient(elem)*elem;
        end
        basis = basis + fTmp*x[i];
    end
    basis = monomials(basis);

    return basis

end


"""
Return the identity matrix of the same dimension and Type as x.
"""
function eye(x)
    n = length(x);
    v = [1.0+0.0x[1]];
    for i = 1:n-1
        push!(v,1.0+0.0x[1]);
    end
    return diagm(v)
end


"""
Plot the first two dimensions of the landscape U. A quiver plot of the vector
field f may optionally be included.
"""
function plotlandscape(f, U, x, lims, vectors=false, scl=0.05)

    Ng = 100;
    Nds = 10;

    # Assemble the grid of points
    xv = collect(linspace(lims[1][1],lims[1][2],Ng));
    yv = collect(linspace(lims[2][1],lims[2][2],Ng));

    # Evaluate U using an array comprehension then plot
    for ii=1:length(x)-2
        U = subs(U,x[ii+2]=>0.0);
        f = subs(f,x[ii+2]=>0.0);
    end
    Umat = [Float64(subs(U, x[1]=>xv[ii], x[2]=>yv[jj])) for ii=1:Ng, jj=1:Ng];
    plt = Plots.contour(xv,yv,Umat'-minimum(Umat),
                        xlabel="x1",ylabel="x2", aspect_ratio=:equal);

    # If desired, evaluate f using an array comprehension then plot
    gU = differentiate(U,x); # Gradient
    fU = f + gU;             # Curl component
    if vectors
        xm = vec([xv[ii] for jj=1:Nds:Ng, ii=1:Nds:Ng]);
        ym = vec([yv[ii] for ii=1:Nds:Ng, jj=1:Nds:Ng]);
        fMat = vec([(scl.*Float64(subs(f[1], x[1]=>xv[ii], x[2]=>yv[jj])),
                     scl.*Float64(subs(f[2], x[1]=>xv[ii], x[2]=>yv[jj])))
               for jj=1:Nds:Ng, ii=1:Nds:Ng]);
        gMat = vec([(scl.*Float64(subs(-gU[1], x[1]=>xv[ii], x[2]=>yv[jj])),
                     scl.*Float64(subs(-gU[2], x[1]=>xv[ii], x[2]=>yv[jj])))
               for jj=1:Nds:Ng, ii=1:Nds:Ng]);
        cMat = vec([(scl.*Float64(subs(fU[1], x[1]=>xv[ii], x[2]=>yv[jj])),
                     scl.*Float64(subs(fU[2], x[1]=>xv[ii], x[2]=>yv[jj])))
               for jj=1:Nds:Ng, ii=1:Nds:Ng]);
        Plots.quiver!(xm,ym, quiver=gMat, color=:black);
        Plots.quiver!(xm,ym, quiver=cMat, color=:black);
    end

    return plt

end


"""
Generate a quiver plot of the vector field f.
"""
function plotvectors(f, x, lims, scl=0.05)

    Ng = 100;
    Nds = 10;

    # Assemble the grid of points
    xv = collect(linspace(lims[1][1],lims[1][2],Ng));
    yv = collect(linspace(lims[2][1],lims[2][2],Ng));

    xm = vec([xv[ii] for jj=1:Nds:Ng, ii=1:Nds:Ng]);
    ym = vec([yv[ii] for ii=1:Nds:Ng, jj=1:Nds:Ng]);
    fMat = vec([(scl.*Float64(subs(f[1], x[1]=>xv[ii], x[2]=>yv[jj])),
                 scl.*Float64(subs(f[2], x[1]=>xv[ii], x[2]=>yv[jj])))
           for jj=1:Nds:Ng, ii=1:Nds:Ng]);
    plt = Plots.quiver(xm,ym, quiver=fMat, color=:black, aspect_ratio=:equal);

    return plt

end


"""
Generate a metric that quantifies the orthogonality of fᵤ with ∇U.
"""
function checknorm(f, U, x)

    # Evaluate ∇U
    ∇U = differentiate(U,x);
    res = dot(∇U,f) + dot(∇U,∇U);

    return mean(abs.(coefficients(res))) /
            mean([mean(abs.(coefficients(f[ii]))) for ii=1:length(x)]);
end

end
