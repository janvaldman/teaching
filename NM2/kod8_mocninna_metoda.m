n=5;                                        % velikost matice
h=1/(n+1); h=1;                              % parametr diskretizace
A=(1/h^2)*gallery('tridiag',n,-1,2,-1);      % sestaveni ridke matice
u=zeros(n,1); u(1)=1;                        % pocatecni vektor
inverse=true;                               % false: mocninna metoda
% true: inverzni metoda
iterMax=1e5;                                 % maximalni pocet iteraci
tol=1e-6;                                    % tolerance
for iter=1:iterMax
    if inverse
        y=A\u;                               % inverzni mocninna metoda
    else
        y=A*u;                               % mocninna metoda
    end
    if norm(y,2)==0
        error('Nulovy vektor: iteraci nelze provest.');
    end
    u=y/norm(y,2);                            % normalizace vektoru
    lambda=u'*(A*u);                          % Rayleighuv podil
    r=norm(A*u-lambda*u,2);                   % reziduum
    fprintf('iter=%d, lambda=%.8f, residual=%.3e\n', ...
        iter,lambda,r);
    fprintf('   u = (');
    fprintf('%.6f ',u);
    fprintf(')^T\n');
    if r<tol
        break
    end
end
if r>=tol
    warning('Dosazen maximalni pocet iteraci bez konvergence.');
end