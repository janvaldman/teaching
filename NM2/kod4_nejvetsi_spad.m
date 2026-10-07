n=4; iters=15;                      % velikost matice a pocet iteraci
A=(n+1)^2*gallery('tridiag',n,-1,2,-1);
b=ones(n,1);
E=@(x) 0.5*dot(x,A*x)-dot(b,x);     % energie
x_SD=zeros(n,iters+1);              % aproximace reseni
r=b-A*x_SD(:,1);                    % reziduum
fprintf('k = %2d, ||r|| = %.4e, E = %.12e\n', ...
    0,norm(r),E(x_SD(:,1)));
for k=1:iters
    Ar=A*r;
    alpha=dot(r,r)/dot(r,Ar);        % delka kroku
    x_SD(:,k+1)=x_SD(:,k)+alpha*r;
    r=r-alpha*Ar;                   % aktualizace rezidua
    fprintf('k = %2d, ||r|| = %.4e, E = %.12e\n', ...
        k,norm(r),E(x_SD(:,k+1)));
end