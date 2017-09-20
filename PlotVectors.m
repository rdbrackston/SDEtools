function [ FigHand ] = PlotVectors( f, U, vars, xlims, ylims )
%PlotLandscape Function to plot the obtained landscape in one or
%two-dimensions
%   Detailed explanation goes here

% Initialise some variables
n = length(vars);

% Evaluate grad(U) then the remaining component
for iv=1:length(vars)
    gU(iv) = diff(U,vars(iv));
end
fU = f + gU';

% Decide on the axis limits and generate grid points
xN = 20;    yN = 20;
[X,Y] = meshgrid(linspace(xlims(1),xlims(2),xN),...
                 linspace(ylims(1),ylims(2),yN));

% Assemble the array over the first two dimensions
Umat = zeros(size(X));
gUx = zeros(size(X));    gUy = zeros(size(X));
fUx = zeros(size(X));    fUy = zeros(size(X));
for ii=1:xN
    for jj=1:yN
        Umat(jj,ii) = subs(U, vars, [X(jj,ii);Y(jj,ii);zeros(n-2,1)]);
        gUx(jj,ii) = subs(gU(1), vars, [X(jj,ii);Y(jj,ii);zeros(n-2,1)]);
        gUy(jj,ii) = subs(gU(2), vars, [X(jj,ii);Y(jj,ii);zeros(n-2,1)]);
        fUx(jj,ii) = subs(fU(1), vars, [X(jj,ii);Y(jj,ii);zeros(n-2,1)]);
        fUy(jj,ii) = subs(fU(2), vars, [X(jj,ii);Y(jj,ii);zeros(n-2,1)]);
    end
end

% Plot the figure
FigHand = figure();
pcolor(X,Y,Umat-double(subs(U,vars,zeros(n,1))))
hold on
quiver([X;X],[Y;Y],[-gUx;fUx],[-gUy;fUy], 'k')
% quiver(X,Y,fUx,fUy, 'w')
set(gca,'TickLabelInterpreter','Latex', 'FontSize',10)
xlabel('$x_1$', 'FontSize',14, 'Interpreter','Latex')
ylabel('$x_2$', 'FontSize',14, 'Interpreter','Latex')
zlabel('$U$', 'FontSize',14, 'Interpreter','Latex')
axis equal

end

