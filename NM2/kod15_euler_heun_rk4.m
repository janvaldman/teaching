%% Explicitni Eulerova, Heunova a Rungeova-Kuttova metoda (RK4)
% Modelova uloha:  y(x)*y'(x) = 2*x^4,  y(0) = 1,  x in [0, 2]
% tj. y' = f(x,y) = 2*x^4/y, presne reseni y(x) = sqrt(4*x^5/5 + 1).
% Pro ruzne pocty kroku n vypiseme chybu v koncovem bode,
% odhadneme rad konvergence a chyby vykreslime v log-log grafu.
clear; close all; clc;

f      = @(x,y) 2*x.^4./y;            % prava strana
yExact = @(x) sqrt(4*x.^5/5 + 1);     % presne reseni
x0 = 0; y0 = 1; xMax = 2;             % pocatecni podminka a konec intervalu
allN = [10 100 1000 10000];           % pocty kroku

metody = {@euler, @heun, @rk4};
nazvy  = {'explicitni Euler', 'Heun', 'RK4'};
chyby  = zeros(numel(metody), numel(allN));

for m = 1:numel(metody)
    fprintf('\n%s\n', nazvy{m});
    for j = 1:numel(allN)
        [x, y] = metody{m}(f, x0, xMax, y0, allN(j));
        h = x(2) - x(1);
        chyby(m,j) = abs(y(end) - yExact(xMax));
        fprintf('   h=%8.2e,  y_h(%g)=%.6f,  e=%10.4e', h, xMax, y(end), chyby(m,j));
        if j > 1                       % experimentalni rad konvergence
            p = log(chyby(m,j-1)/chyby(m,j)) / log(allN(j)/allN(j-1));
            fprintf(',  p=%5.2f', p);
        end
        fprintf('\n');
    end
end

% Pozn.: u RK4 je chyba pro n = 10000 uz na urovni zaokrouhlovacich chyb
% (presnost typu double), proto tam odhad radu p ztraci smysl.

% Log-log graf chyb (body na urovni zaokrouhlovacich chyb vynechame)
h = (xMax - x0)./allN;
figure;
for m = 1:numel(metody)
    ok = chyby(m,:) > 1e-13;
    loglog(h(ok), chyby(m,ok), '-o', 'LineWidth', 1.5); hold on;
end
hold off; grid on; set(gca, 'XDir', 'reverse');
xlabel('h'); ylabel('e(x_{max})');
legend('explicitni Euler (p=1)', 'Heun (p=2)', 'RK4 (p=4)', 'Location', 'southwest');
title('Chyba v koncovem bode');

% Aproximace pro n = 10 a presne reseni
figure;
[xE, yE] = euler(f, x0, xMax, y0, 10);
[xH, yH] = heun(f, x0, xMax, y0, 10);
[xR, yR] = rk4(f, x0, xMax, y0, 10);
xx = linspace(x0, xMax, 200);
plot(xx, yExact(xx), 'k-', xE, yE, 'o-', xH, yH, 's-', xR, yR, '^-');
legend('presne reseni', 'explicitni Euler', 'Heun', 'RK4', 'Location', 'northwest');
xlabel('x'); ylabel('y'); title('Aproximace pro n = 10'); grid on;

%% Jednotlive metody
function [x, y] = euler(f, x0, xMax, y0, n)      % explicitni Eulerova metoda
    x = linspace(x0, xMax, n+1)'; h = x(2) - x(1);
    y = zeros(n+1, 1); y(1) = y0;
    for i = 1:n
        y(i+1) = y(i) + h*f(x(i), y(i));
    end
end

function [x, y] = heun(f, x0, xMax, y0, n)       % Heunova metoda
    x = linspace(x0, xMax, n+1)'; h = x(2) - x(1);
    y = zeros(n+1, 1); y(1) = y0;
    for i = 1:n
        k1 = f(x(i), y(i));
        k2 = f(x(i+1), y(i) + h*k1);
        y(i+1) = y(i) + h/2*(k1 + k2);
    end
end

function [x, y] = rk4(f, x0, xMax, y0, n)        % Rungeova-Kuttova metoda 4. radu
    x = linspace(x0, xMax, n+1)'; h = x(2) - x(1);
    y = zeros(n+1, 1); y(1) = y0;
    for i = 1:n
        k1 = f(x(i),       y(i));
        k2 = f(x(i) + h/2, y(i) + h/2*k1);
        k3 = f(x(i) + h/2, y(i) + h/2*k2);
        k4 = f(x(i) + h,   y(i) + h*k3);
        y(i+1) = y(i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end
