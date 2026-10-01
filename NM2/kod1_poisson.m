close all

n = 20; h = 1/(n + 1);   % Parametry site
A = (1/h^2)*gallery('tridiag',n,-1,2,-1);  % Sestaveni ridke matice
x = h*linspace(1,n,n)';    % Uzly uvnitr intervalu
b = x.^2;   % Prava strana

u = A \ b;  % Numericke reseni


uExact = @(y) (y-y.^4)/12;  % Presne reseni

% Graf numerickeho a presneho reseni
scatter(x,u,'b','MarkerFaceColor','b')
hold on
xx = linspace(0,1);
plot(xx,uExact(xx),'r-')
hold off

error = max(abs(u-uExact(x)))   % Maximalni chyba