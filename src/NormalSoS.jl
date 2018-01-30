module NormalSoS

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials, CSDP
using Plots
gr()
export normdecomp

# To do:
#  - Add plotting functions
#  - Add function to check normality

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


function normopt2(f, x, basis, SDPsolver=CSDPSolver())
    # Vector ϵ used for lower bound

    n = length(f);

    m = SOSModel(solver=SDPsolver);

    @polyvariable m V basis

    # Specify the lower bounding polynomial, bnd
    o = zeros(Int,1,n)
    o[1] = maximum([degree(Vi,x[1]) for Vi in basis])
    bnd = x[1]^o[1];
    if n>1
        # First add the maximum even order polynomial for each x[ii]
        for ii=2:n;
            o[ii] = maximum([degree(Vi,x[ii]) for Vi in basis])
            if o[ii]%2==0
                bnd += x[ii]^o[ii];
            end
        end
        # Now add any fully even terms of mixed x[ii], e.g. x[1]^2*x[2]^2
        for bi in basis
            e = exponents(bi);
            if all(y->y==e[1]&&e[1]%2==0, e) && e[1]!=0
                bnd += bi
            end
        end
        b = length(bnd);
    else
        b = 1;
        bnd = [bnd];
    end

    # Positive definitness constraint on V
    @variable m ϵ[1:b];
    for ii=1:b; @constraint m ϵ[ii] ≥ 0; end
    # @constraint m sum([ϵ[ii]*bnd[ii] for ii=1:length(bnd)]) ≥ 0
    @polyconstraint m V ≥ sum([ϵ[ii]*bnd[ii] for ii=1:b])

    # Apply matrix constraint, ∇U⋅g ≤ 0.
    I = NormalSoS.eye(x);
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
    @SDconstraint m Mv ⪰ 0 # Mv positive definite

    @objective m Max sum(ϵ)
    status = solve(m);

    @show status
    @show bnd
    @show getvalue(ϵ)
    return getvalue(V)

end


function lyapunov(f, x, SDPsolver=CSDPSolver(), o=2)

    n = length(f);

    m = SOSModel(solver=SDPsolver);
    @variable m ϵ
    @polyvariable m V monomials(x,o);

    # Positive definiteness constraint
    @polyconstraint m V ≥ ϵ*sum(x.^2);
    @constraint m ϵ ≥ 0
    @objective m Max ϵ

    # Standard Lyapunov constraint
    @constraint m -dot(differentiate(V, x),f) ≥ 0;

    status = solve(m)
    @show status
    @show getvalue(ϵ)
    return getvalue(V)

end


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

function eye(x)
    n = length(x);
    v = [1.0+0.0x[1]];
    for i = 1:n-1
        push!(v,1.0+0.0x[1]);
    end
    return diagm(v)
end

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

function checknorm(f, U, x)

    # Evaluate ∇U
    ∇U = differentiate(U,x);
    res = dot(∇U,f) + dot(∇U,∇U);

    return mean(abs.(coefficients(res))) /
            mean([mean(abs.(coefficients(f[ii]))) for ii=1:length(x)]);
end

end
