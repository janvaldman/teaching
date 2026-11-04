clear;
clc;
close all;

n = 10;
A = full(gallery('tridiag', n, -1, 2, -1));
A(1,2) = -3;                    % matice je nyní nesymetrická

tol = 1e-4;                     % tolerance pro poddiagonální prvky
maxIter = 1000;                 % maximální počet iterací
Asim = A;
V = eye(n);
iter = 0;

while norm(tril(Asim, -1), 'fro') > tol && iter < maxIter
    iter = iter + 1;

    mu = 0*Asim(n,n);             % jednoduchý posun z pravého dolního prvku
    [Q, R] = qr(Asim - mu*eye(n));
    Asim = R*Q + mu*eye(n);
    V = V*Q;
end

chybaPodDiagonaly = norm(tril(Asim, -1), 'fro');
fprintf('Počet iterací: %d\n', iter);

if chybaPodDiagonaly > tol
    warning('QR iterace nedosáhla tolerance v maxIter krocích.');
end

% Aproximace Schurovy horní trojúhelníkové matice
T = triu(Asim);

% Aproximace vlastních čísel z diagonály
vlastniCislaQR = diag(T)

% Kontrola vztahu A*V = V*Asim
kontrolaPodobnosti = norm(A*V - V*Asim, 'fro')

% Vlastní vektory T a jejich převod zpět k původní matici A
[Y, D] = eig(T);
W = V*Y;

% Kontrola vlastních vektorů původní matice
kontrolaVlastnichVektoru = norm(A*W - W*D, 'fro')