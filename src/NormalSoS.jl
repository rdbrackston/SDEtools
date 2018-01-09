module NormalSoS

using SumOfSquares, JuMP, PolyJuMP, DynamicPolynomials, MultivariatePolynomials, CSDP
using Plots
gr()
export normdecomp

# To do:
#  - Add plotting functions
#  - Add function to check normality

function normdecomp(f, x, SDPsolver=CSDPSolver(), nIters=1, o=2, basis=:minimal)

    n = length(f);

    m1 = SOSModel(solver=SDPsolver);
    @variable m1 ϵ

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
    @polyvariable m1 V Z
    # @variable m1 V Z

    # Positive definiteness constraint
    @polyconstraint m1 V ≥ ϵ*sum(x.^o);
    # @constraint m1 V ≥ ϵ*sum(x.^o);

    # Apply matrix constraint, ∇U⋅g ≤ 0.
    I = NormalSoS.eye(x);
    Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
    # Mv = [-dot(differentiate(V, x),f) differentiate(V,x)';
    #        differentiate(V,x)         I];
    for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
    @SDconstraint m1 Mv ⪰ 0 # Mv positive definite

    @objective m1 Max ϵ

    status = solve(m1);

    # Now iterate to improve
    for ii=1:nIters

        U = getvalue(V);

        m2 = SOSModel(solver=SDPsolver);
        @variable m2 ϵ
        @variable m2 α
        @polyvariable m2 V Z
        # @variable m2 V Z

        # Positive definiteness constraint
        @polyconstraint m2 V ≥ ϵ*sum(x.^o);
        # @constraint m2 V ≥ ϵ*sum(x.^o);

        # Apply matrix constraint, ∇U⋅g ≤ 0.
        Mv = [-dot(differentiate(V, x),f); differentiate(V,x)];
        # Mv = [-dot(differentiate(V, x),f) differentiate(V,x)';
        #        differentiate(V,x)         I];
        for ii=1:n; Mv = hcat(Mv, [differentiate(V,x)[ii];I[:,ii]]); end
        @SDconstraint m2 Mv ⪰ 0 # Mv positive definite

        # Wynn inequality constraint
        @polyconstraint m2 dot(differentiate(V,x),f+2*differentiate(U,x)) ≥
            α*dot(differentiate(U,x),f) + (1+α)*sum(differentiate(U,x).^2)
        # @constraint m2 dot(differentiate(V,x),f+2*differentiate(U,x)) ≥
        #     α*dot(differentiate(U,x),f) + (1+α)*sum(differentiate(U,x).^2)

        @constraint m2 α ≥ 0
        @constraint m2 ϵ ≥ 0
        @objective m2 Min α

        status = solve(m2);

    end

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
    plt = Plots.contour(xv,yv,Umat', xlabel="x1",ylabel="x2", aspect_ratio=:equal);

    # If desired, evaluate f using an array comprehension then plot
    if vectors
        xm = vec([xv[ii] for jj=1:Nds:Ng, ii=1:Nds:Ng]);
        ym = vec([yv[ii] for ii=1:Nds:Ng, jj=1:Nds:Ng]);
        fMat = vec([(scl.*Float64(subs(f[1], x[1]=>xv[ii], x[2]=>yv[jj])),
                     scl.*Float64(subs(f[2], x[1]=>xv[ii], x[2]=>yv[jj])))
               for jj=1:Nds:Ng, ii=1:Nds:Ng]);
        Plots.quiver!(xm,ym, quiver=fMat, color=:black);
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
