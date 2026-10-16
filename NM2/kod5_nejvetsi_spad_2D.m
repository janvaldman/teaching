clear; clc; close all;

A = [4 1; 1 2];
b = [1; 1];
iters = 15;
x0 = [0; 0];

E = @(x) 0.5*dot(x,A*x) - dot(b,x);

% Metoda nejvetsiho spadu
x_SD = zeros(2,iters+1);
x_SD(:,1) = x0;
r = b-A*x0;

fprintf(' k       x_1          x_2          ||r||          E(x)\n');
fprintf('---------------------------------------------------------\n');
fprintf('%2d   %10.6f   %10.6f   %.4e   %.8f\n', ...
    0,x0(1),x0(2),norm(r),E(x0));

for k = 1:iters
    Ar = A*r;
    alpha = dot(r,r)/dot(r,Ar);

    x_SD(:,k+1) = x_SD(:,k) + alpha*r;
    r = r-alpha*Ar;

    fprintf('%2d   %10.6f   %10.6f   %.4e   %.8f\n', ...
        k,x_SD(1,k+1),x_SD(2,k+1),norm(r),E(x_SD(:,k+1)));

    if norm(r) < 1e-12
        x_SD = x_SD(:,1:k+1);
        break
    end
end

% Presne reseni
x = A\b;
fprintf('\nPresne reseni: x = (%.6f, %.6f)^T\n',x(1),x(2));
fprintf('Minimalni energie: E(x) = %.8f\n',E(x));

% Automaticke urceni mezi grafu
X = [x_SD, x];

xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));

margin = 0.2*max([xmax-xmin, ymax-ymin, 0.1]);

xmin = xmin-margin;
xmax = xmax+margin;
ymin = ymin-margin;
ymax = ymax+margin;

% Mrizka pro vykresleni vrstevnic
[x1,x2] = meshgrid(linspace(xmin,xmax,200), ...
                   linspace(ymin,ymax,200));

% Energeticka funkce
F = 0.5*(A(1,1)*x1.^2 + 2*A(1,2)*x1.*x2 ...
    + A(2,2)*x2.^2) - b(1)*x1 - b(2)*x2;

% Graf vrstevnic
figure;
contour(x1,x2,F,20);
hold on;

% Trajektorie iteraci
plot(x_SD(1,:),x_SD(2,:),'ro-', ...
    'LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',5);

% Pocatecni aproximace a presne reseni
plot(x0(1),x0(2),'bo', ...
    'MarkerFaceColor','b','MarkerSize',8);

plot(x(1),x(2),'p', ...
    'MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',12);

% Popisky iteraci
for k = 0:size(x_SD,2)-1
    text(x_SD(1,k+1)+0.02*margin, ...
         x_SD(2,k+1)+0.02*margin, ...
         sprintf('x_%d',k),'FontSize',9);
end

xlabel('x_1');
ylabel('x_2');
title('Metoda nejvetsiho spadu');

legend('Vrstevnice E(x)','Iterace','x_0','Presne reseni', ...
    'Location','southeast');

axis equal;
xlim([xmin xmax]);
ylim([ymin ymax]);
grid on;
hold off;

% Ulozeni obrazku do slozky skriptu
scriptFolder = fileparts(mfilename('fullpath'));

if isempty(scriptFolder)
    scriptFolder = pwd;
end

filename = fullfile(scriptFolder,'steepest_descent.pdf');
exportgraphics(gcf,filename,'ContentType','vector');

fprintf('\nObrazek ulozen: %s\n',filename);