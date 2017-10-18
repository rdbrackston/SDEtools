function [ FigHand ] = PlotLandscape( f, U, vars, xlims, ylims )
%PlotLandscape Function to plot the obtained landscape in one or
%two-dimensions
%   Detailed explanation goes here

% Initialise some variables
n = length(vars);

% Decide on the axis limits and generate grid points
xN = 30;    yN = 30;
[Y,X] = meshgrid(linspace(xlims(1),xlims(2),xN),...
                 linspace(ylims(1),ylims(2),yN));

% Assemble the array over the first two dimensions
Umat = zeros(size(X));
for ii=1:xN
    for jj=1:yN
        Umat(jj,ii) = subs(U, vars, [X(jj,ii);Y(jj,ii);zeros(n-2,1)]);
    end
end

% Plot the figure
FigHand = figure();
surf(X,Y,Umat-double(subs(U,vars,zeros(n,1))))
set(gca,'TickLabelInterpreter','Latex', 'FontSize',10)
xlabel('$x_1$', 'FontSize',14, 'Interpreter','Latex')
ylabel('$x_2$', 'FontSize',14, 'Interpreter','Latex')
zlabel('$U$', 'FontSize',14, 'Interpreter','Latex')

end

