%============================================================
% Mocninne metody pro vypocet vlastniho cisla
%============================================================
clc
close all

% Matice
A = [4 1; 2 3];                 % matice z prednasek 5 a 6
u = [0; 1];                     % pocatecni vektor

% Volba metody
inverse = true;                 % false: mocninna, true: inverzni
shift = 2.2;                    % 0: bez posunu, jinak posun sigma

eigenvalueMethod = 'components'; % 'components' nebo 'Rayleigh'

iterMax = 50;                   % maximalni pocet iteraci
tol = 1e-6;                     % tolerance rezidua

% Posunuta matice
B = A - shift*eye(size(A));

if inverse || shift ~= 0
    if rcond(full(B)) < eps
        error('Matice A-shift*I je singularni nebo spatne podminena.');
    end
end

for iter = 1:iterMax

    % Iteracni krok
    if inverse
        y = B\u;                % inverzni nebo posunuta inverzni metoda
    else
        y = B*u;                % mocninna nebo posunuta mocninna metoda
    end

    if norm(y,2) == 0
        error('Nulovy vektor: iteraci nelze provest.');
    end

    u = y/norm(y,2);            % normalizace vektoru

    % Soucin pro odhad vlastniho cisla a reziduum
    Au = A*u;

    switch eigenvalueMethod
        case 'components'
            [~,i] = max(abs(u));
            lambda = Au(i)/u(i);

        case 'Rayleigh'
            lambda = u'*Au;     % u je normalizovany

        otherwise
            error('Neznamy zpusob vypoctu vlastniho cisla.');
    end

    % Reziduum vlastni dvojice pro puvodni matici A
    r = norm(Au - lambda*u,2);

    fprintf('iter=%d, lambda=%.8f, residual=%.3e\n', ...
        iter, lambda, r);
    fprintf('   u = (');
    fprintf('%.6f ', u);
    fprintf(')^T\n');

    if r < tol
        break
    end
end

if r >= tol
    warning('Dosazen maximalni pocet iteraci bez konvergence.');
end
