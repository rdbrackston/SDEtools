function [ metric ] = CheckNorm( f, U, vars )
%CheckNorm Function to check the normality of the decomposition
%   Detailed explanation goes here

% Evaluate grad(U) and coefficients in f
for iv=1:length(vars)
    gU(iv) = diff(U,vars(iv));
    minCoF(iv) = min(abs(coeffs(f(iv))));
end

% Evaluate gradU \cdot gU
res = gU*f + gU*gU.';

% Check coefficients relative to those in f

metric = max(abs(coeffs(res))) / min(minCoF);

end

