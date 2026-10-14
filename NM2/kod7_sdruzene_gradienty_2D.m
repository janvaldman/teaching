clear; clc; close all;

n = 2;                              % rozmer soustavy
iters = 15;                         % maximalni pocet iteraci
tol = 1e-12;                        % tolerance

% Matice a prava strana
A = [4 1; 1 2];
b = ones(n,1);

E = @(x) 0.5*dot(x,A*x)-dot(b,x);  % energie

% Metoda sdruzenych gradientu
x_CG = zeros(n,iters+1);           % aproximace reseni
r = b-A*x_CG(:,1);                 % reziduum
p = r;                             % pocatecni smer

fprintf(' k       ||r||          E(x)\n');
fprintf('----------------------------------\n');
fprintf('%2d   %.4e   %.12e\n',0,norm(r),E(x_CG(:,1)));

for k = 1:iters
    Ap = A*p;
    alpha = dot(r,r)/dot(p,Ap);

    x_CG(:,k+1) = x_CG(:,k)+alpha*p;
    r_new = r-alpha*Ap;

    fprintf('%2d   %.4e   %.12e\n', ...
        k,norm(r_new),E(x_CG(:,k+1)));

    if norm(r_new) < tol
        x_CG = x_CG(:,1:k+1);
        break
    end

    beta = dot(r_new,r_new)/dot(r,r);
    p = r_new+beta*p;
    r = r_new;
end

% Presne reseni
x = A\b;
fprintf('\nPresne reseni: x = (%.6f, %.6f)^T\n',x(1),x(2));

% Vizualizace pro n = 2
if n == 2
    [x1,x2] = meshgrid(linspace(-0.1,0.6,200), ...
        linspace(-0.1,0.6,200));

    F = 0.5*(A(1,1)*x1.^2 + 2*A(1,2)*x1.*x2 ...
        + A(2,2)*x2.^2) - b(1)*x1 - b(2)*x2;

    figure;
    contour(x1,x2,F,20);
    hold on;

    plot(x_CG(1,:),x_CG(2,:),'ro-', ...
        'LineWidth',1.5,'MarkerFaceColor','r');

    plot(x_CG(1,1),x_CG(2,1),'bo', ...
        'MarkerFaceColor','b','MarkerSize',8);

    plot(x(1),x(2),'p','MarkerFaceColor','g', ...
        'MarkerEdgeColor','g','MarkerSize',12);

    for k = 0:size(x_CG,2)-1
        text(x_CG(1,k+1)+0.008,x_CG(2,k+1)+0.008, ...
            sprintf('x_%d',k));
    end

    xlabel('x_1');
    ylabel('x_2');
    title('Metoda sdruzenych gradientu (CG)');
    legend('Vrstevnice E(x)','Iterace CG','x_0', ...
        'Presne reseni','Location','best');
    axis equal;
    grid on;
    hold off;
end