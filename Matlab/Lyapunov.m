function [ U ] = Lyapunov( f, vars, o )
%Lyapunov Function to find a Lyapunov function via the standard method
%   This function tries to find a Lyapunov function U(x) for the vector
%   field f(x) using the method oulined in the SOSTOOLS manual. Such a
%   Lyapunov function is generally not unique and will not be orthogonal
%   but may provide a useful starting point in some cases.

syms epsl

prog = sosprogram(vars,epsl);
[prog,V] = sospolyvar(prog,vars.^o,'wscoeff');

% =============================================
% Constraint 1 : positive definiteness of V
% V(x) - epsilon(x1^2 + x2^2 + ... xn^2) >= 0
p = 2;
prog = sosineq(prog,V-epsl*sum(vars.^p));

% Constraint 2: decreasing V along systrem trajectories
% grad(V).f <= 0
for iv=1:length(vars)
    gV(iv) = diff(V,vars(iv));
end
expr = -gV*f;
prog = sosineq(prog,expr);

% Constraint 3: epsilon positive
prog = sosineq(prog,epsl)

% =============================================
% Set objective to maximise epsilon then call solver
prog = sossetobj(prog, -epsl);
prog = sossolve(prog);

U = sosgetsol(prog,V);

end

