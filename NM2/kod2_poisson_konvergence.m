close all

% Pocty vnitrnich uzlu
nValues = 10*2.^(0:9)-1;

% Presne reseni
uExact = @(y) (y-y.^4)/12;

errors = zeros(size(nValues));
hValues = zeros(size(nValues));

for k = 1:length(nValues)

    % Parametry site
    n = nValues(k);
    h = 1/(n + 1);

    % Sestaveni ridke matice
    A = (1/h^2)*gallery('tridiag',n,-1,2,-1);

    % Uzly uvnitr intervalu
    x = h*linspace(1,n,n)';

    % Prava strana
    b = x.^2;

    % Numericke reseni
    u = A \ b;

    % Maximalni chyba
    errors(k) = max(abs(u-uExact(x)));
    hValues(k) = h;

end

% Vypis vysledku
disp(table(hValues',errors', ...
    'VariableNames',{'h','error'}))

% Graf chyby
loglog(hValues,errors,'o-')
set(gca,'XDir','reverse')
grid on
xlabel('h')
ylabel('maximalni chyba')
title('Konvergence diferencni metody')

% Experimentalni rad konvergence
p = log(errors(1:end-1)./errors(2:end))/log(2);

% Vypis experimentalniho radu
disp(table(hValues(2:end)',p', ...
    'VariableNames',{'h','order'}))