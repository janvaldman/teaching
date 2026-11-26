% Explicitni a implicitni Euler pro ulohu y*y' = 2*x^4, y(0) = 1 na [0,2]
f = @(x,y) 2*x.^4./y;                   % prava strana y' = f(x,y)
dfdy = @(x,y) -2*x.^4./y.^2;            % derivace f podle y (pro Newtona)
yExact = @(x) sqrt(4*x.^5/5 + 1);       % presne reseni
x0 = 0; y0 = 1; xMax = 2;               % pocatecni podminka a konec intervalu
nAll = [10 100 1000 10000];             % pocty kroku

e = zeros(numel(nAll),2);               % chyby v bode xMax: explicitni, implicitni
for j = 1:numel(nAll)
    n = nAll(j); h = (xMax - x0)/n; x = x0 + h*(0:n);
    yE = y0*ones(1,n+1); yI = yE;       % aproximace v uzlech x
    for i = 1:n
        yE(i+1) = yE(i) + h*f(x(i),yE(i));               % explicitni Euler

        z = yI(i);                                       % implicitni Euler:
        for k = 1:20                                     % Newton pro
            g  = z - yI(i) - h*f(x(i+1),z);              % g(z) = z - y_i - h*f(x_{i+1},z) = 0
            dg = 1 - h*dfdy(x(i+1),z);
            z  = z - g/dg;
        end
        yI(i+1) = z;
    end
    e(j,:) = abs([yE(end) yI(end)] - yExact(xMax));

    if j == 1                           % graf reseni pro nejmensi n
        figure(1)
        xx = linspace(x0, xMax);
        plot(xx, yExact(xx), 'k-', x, yE, 'o-', x, yI, 's-')
        grid on; xlabel('x'); ylabel('y')
        legend('presne reseni', 'explicitni Euler', 'implicitni Euler', 'Location', 'northwest')
        title(sprintf('Aproximace reseni pro n = %d', n))
    end
end

% Rad konvergence: krok se zmensuje 10x, proto p = log10(e_h/e_{h/10})
p = log10(e(1:end-1,:)./e(2:end,:));

% Vypis tabulky
hAll = (xMax - x0)./nAll;
fprintf('\n  %7s |  %9s  %5s |  %9s  %5s\n', 'h', 'explicitni', 'p', 'implicitni', 'p');
fprintf('%s\n', repmat('-', 1, 46));
for j = 1:numel(nAll)
    fprintf('  %7.1e', hAll(j));
    for m = 1:2
        if j == 1, pStr = '-'; else, pStr = sprintf('%.2f', p(j-1,m)); end
        fprintf(' |  %9.2e  %5s', e(j,m), pStr);
    end
    fprintf('\n');
end
fprintf('(explicitni metoda reseni podhodnocuje, implicitni nadhodnocuje)\n\n');
