% Explicitni Euler, Heun a RK4 pro ulohu y*y' = 2*x^4, y(0) = 1 na [0,2]
f = @(x,y) 2*x.^4./y;                   % prava strana y' = f(x,y)
yExact = @(x) sqrt(4*x.^5/5 + 1);       % presne reseni
x0 = 0; y0 = 1; xMax = 2;               % pocatecni podminka a konec intervalu
nAll = [10 100 1000 10000];             % pocty kroku

e = zeros(numel(nAll),3);              % chyby v bode xMax: Euler, Heun, RK4
for j = 1:numel(nAll)
    n = nAll(j); h = (xMax - x0)/n; x = x0 + h*(0:n);
    yE = y0*ones(1,n+1); yH = yE; yR = yE;  % aproximace v uzlech x
    for i = 1:n
        yE(i+1) = yE(i) + h*f(x(i),yE(i));                    % explicitni Euler

        k1 = f(x(i),yH(i)); k2 = f(x(i+1),yH(i) + h*k1);      % Heun
        yH(i+1) = yH(i) + h/2*(k1 + k2);

        k1 = f(x(i),yR(i));           k2 = f(x(i)+h/2,yR(i) + h/2*k1);   % RK4
        k3 = f(x(i)+h/2,yR(i) + h/2*k2); k4 = f(x(i+1),yR(i) + h*k3);
        yR(i+1) = yR(i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    e(j,:) = abs([yE(end) yH(end) yR(end)] - yExact(xMax));

    xx = linspace(x0, xMax);
    if j == 1                           % obrazek 1: jen explicitni Euler (n = 10, pozdeji n = 100)
        figure(1)
        plot(xx, yExact(xx), 'k-', x, yE, 'o-', 'LineWidth', 1.5); hold on
        grid on; xlabel('x'); ylabel('y')
        title('Explicitni Euler')

        figure(2)                       % obrazek 2: srovnani metod pro n = 10
        plot(xx, yExact(xx), 'k-', x, yE, 'o-', x, yH, 's-', x, yR, '^-', 'LineWidth', 1.5)
        grid on; xlabel('x'); ylabel('y')
        legend('presne reseni', 'Euler', 'Heun', 'RK4', 'Location', 'northwest')
        title(sprintf('Srovnani metod pro n = %d', n))
    end
    if j == 2                           % do obrazku 1 pridame Eulera pro n = 100
        figure(1)
        plot(x, yE, '--', 'LineWidth', 1.5); hold off
        legend('presne reseni', 'Euler, n = 10', 'Euler, n = 100', 'Location', 'northwest')
    end
end

% Rad konvergence: krok se zmensuje 10x, proto p = log10(e_h/e_{h/10})
p = log10(e(1:end-1,:)./e(2:end,:));

% Vypis tabulky
hAll = (xMax - x0)./nAll;
fprintf('\n  %7s |  %9s  %5s |  %9s  %5s |  %9s  %5s\n', 'h', 'Euler', 'p', 'Heun', 'p', 'RK4', 'p');
fprintf('%s\n', repmat('-', 1, 66));
for j = 1:numel(nAll)
    fprintf('  %7.1e', hAll(j));
    for m = 1:3
        if j == 1, pStr = '-'; else, pStr = sprintf('%.2f', p(j-1,m)); end
        fprintf(' |  %9.2e  %5s', e(j,m), pStr);
    end
    fprintf('\n');
end
fprintf('(u RK4 je chyba pro nejmensi h uz na urovni zaokrouhlovacich chyb)\n\n');

% Graf chyb v logaritmickych souradnicich
figure(3)
loglog(hAll, e(:,1), 'o-', hAll, e(:,2), 's-', hAll(1:3), e(1:3,3), '^-', ...
       'LineWidth', 1.5, 'MarkerSize', 7)
set(gca, 'XDir', 'reverse', 'YTick', 10.^(-12:3:0)); grid on
xlim([1e-4 1]); ylim([1e-13 1])
xlabel('h'); ylabel('chyba v x_{max}')
legend('Euler (p=1)', 'Heun (p=2)', 'RK4 (p=4)', 'Location', 'southwest')
title('Chyba v koncovem bode')
