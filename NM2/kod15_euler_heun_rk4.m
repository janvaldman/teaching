% Explicitni Euler, Heun a RK4 pro ulohu y*y' = 2*x^4, y(0) = 1 na [0,2]
close all
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

% Srovnani pri stejne cene: 200 vyhodnoceni f
% (Euler 1 vyhodnoceni/krok -> n = 200, Heun 2 -> n = 100, RK4 4 -> n = 50)
nCena = [200 100 50]; eCena = zeros(1,3);
for m = 1:3
    n = nCena(m); h = (xMax - x0)/n; x = x0 + h*(0:n); y = y0;
    for i = 1:n
        k1 = f(x(i),y);
        if m == 1                                   % Euler
            y = y + h*k1;
        elseif m == 2                               % Heun
            k2 = f(x(i+1),y + h*k1); y = y + h/2*(k1 + k2);
        else                                        % RK4
            k2 = f(x(i)+h/2,y + h/2*k1); k3 = f(x(i)+h/2,y + h/2*k2);
            k4 = f(x(i+1),y + h*k3); y = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
        end
    end
    eCena(m) = abs(y - yExact(xMax));
end
fprintf('  Stejna cena (200 vyhodnoceni f):\n');
fprintf('    Euler  h = %5.3f   chyba = %9.2e\n', (xMax-x0)/nCena(1), eCena(1));
fprintf('    Heun   h = %5.3f   chyba = %9.2e\n', (xMax-x0)/nCena(2), eCena(2));
fprintf('    RK4    h = %5.3f   chyba = %9.2e\n\n', (xMax-x0)/nCena(3), eCena(3));

% Odhad radu BEZ znalosti presneho reseni: tri vypocty s kroky h, h/2, h/4
%   y_h - y_{h/2} ~ C h^p (1 - 2^(-p))  =>  p ~ log2(|y_h - y_{h/2}| / |y_{h/2} - y_{h/4}|)
nOdh = [10 20 40];                      % h = 0.2, 0.1, 0.05
yOdh = zeros(3,3);                      % radky: metody, sloupce: h, h/2, h/4
for m = 1:3
    for k = 1:3
        n = nOdh(k); h = (xMax - x0)/n; x = x0 + h*(0:n); y = y0;
        for i = 1:n
            k1 = f(x(i),y);
            if m == 1                                   % Euler
                y = y + h*k1;
            elseif m == 2                               % Heun
                k2 = f(x(i+1),y + h*k1); y = y + h/2*(k1 + k2);
            else                                        % RK4
                k2 = f(x(i)+h/2,y + h/2*k1); k3 = f(x(i)+h/2,y + h/2*k2);
                k4 = f(x(i+1),y + h*k3); y = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
            end
        end
        yOdh(m,k) = y;
    end
end
pOdh = log2(abs(yOdh(:,1) - yOdh(:,2))./abs(yOdh(:,2) - yOdh(:,3)));
fprintf('  Odhad radu bez presneho reseni (h = 0.2, 0.1, 0.05):\n');
fprintf('    Euler p = %.2f,  Heun p = %.2f,  RK4 p = %.2f\n\n', pOdh);

% Graf chyb v logaritmickych souradnicich
figure(3)
loglog(hAll, e(:,1), 'o-', hAll, e(:,2), 's-', hAll(1:3), e(1:3,3), '^-', ...
       'LineWidth', 1.5, 'MarkerSize', 7)
set(gca, 'XDir', 'reverse', 'YTick', 10.^(-12:3:0)); grid on
xlim([1e-4 1]); ylim([1e-13 1])
xlabel('h'); ylabel('chyba v x_{max}')
legend('Euler (p=1)', 'Heun (p=2)', 'RK4 (p=4)', 'Location', 'southwest')
title('Chyba v koncovem bode')
