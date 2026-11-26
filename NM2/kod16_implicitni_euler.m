%% Implicitni Eulerova metoda a srovnani s explicitni
% Modelova uloha:  y' = f(x,y) = 2*x^4/y,  y(0) = 1,  x in [0, 2].
% Implicitni krok  y_{i+1} = y_i + h*f(x_{i+1}, y_{i+1})  resime
% Newtonovou metodou pro g(z) = z - y_i - h*f(x_{i+1}, z) = 0.
% Pro tuto ulohu lze krok spocitat i primo (kladny koren kvadraticke rovnice):
%   y_{i+1} = (y_i + sqrt(y_i^2 + 8*h*x_{i+1}^4))/2.
clear; close all; clc;

f      = @(x,y) 2*x.^4./y;            % prava strana
dfdy   = @(x,y) -2*x.^4./y.^2;        % derivace f podle y (pro Newtona)
yExact = @(x) sqrt(4*x.^5/5 + 1);     % presne reseni
x0 = 0; y0 = 1; xMax = 2;
allN = [10 100 1000 10000];

fprintf('       h     explicitni Euler      implicitni Euler (Newton)   kontrola\n');
for n = allN
    [x, yE] = euler(f, x0, xMax, y0, n);
    [~, yI] = implicitni_euler(f, dfdy, x0, xMax, y0, n);
    yQ = implicitni_euler_vzorec(x0, xMax, y0, n);
    h = x(2) - x(1);
    fprintf('%8.1e   %.6f (e=%8.2e)   %.6f (e=%8.2e)   %8.1e\n', h, ...
        yE(end), abs(yE(end)-yExact(xMax)), yI(end), abs(yI(end)-yExact(xMax)), ...
        abs(yI(end)-yQ));
end
% Explicitni metoda reseni podhodnocuje, implicitni nadhodnocuje,
% obe jsou prvniho radu.

% Graf pro n = 10
n = 10;
[x, yE] = euler(f, x0, xMax, y0, n);
[~, yI] = implicitni_euler(f, dfdy, x0, xMax, y0, n);
xx = linspace(x0, xMax, 200);
plot(xx, yExact(xx), 'k-', x, yE, 'o-', x, yI, 's-', x, (yE+yI)/2, 'x--');
legend('presne reseni', 'explicitni Euler', 'implicitni Euler', 'prumer obou', ...
    'Location', 'northwest');
xlabel('x'); ylabel('y'); title('Explicitni a implicitni Euler, n = 10'); grid on;

%% Jednotlive metody
function [x, y] = euler(f, x0, xMax, y0, n)      % explicitni Eulerova metoda
    x = linspace(x0, xMax, n+1)'; h = x(2) - x(1);
    y = zeros(n+1, 1); y(1) = y0;
    for i = 1:n
        y(i+1) = y(i) + h*f(x(i), y(i));
    end
end

function [x, y] = implicitni_euler(f, dfdy, x0, xMax, y0, n)
    x = linspace(x0, xMax, n+1)'; h = x(2) - x(1);
    y = zeros(n+1, 1); y(1) = y0;
    for i = 1:n
        z = y(i);                                % pocatecni odhad Newtona
        for k = 1:50
            g  = z - y(i) - h*f(x(i+1), z);
            dg = 1 - h*dfdy(x(i+1), z);
            dz = g/dg; z = z - dz;
            if abs(dz) < 1e-14*max(1, abs(z)), break; end
        end
        y(i+1) = z;
    end
end

function y = implicitni_euler_vzorec(x0, xMax, y0, n)   % primy vzorec
    x = linspace(x0, xMax, n+1)'; h = x(2) - x(1);
    y = y0;
    for i = 1:n
        y = (y + sqrt(y^2 + 8*h*x(i+1)^4))/2;
    end
end
