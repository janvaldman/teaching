n=4; iters=15;                      % velikost matice a pocet iteraci
tol=1e-12;                          % tolerance
A=(n+1)^2*gallery('tridiag',n,-1,2,-1);
b=ones(n,1);

E=@(x) 0.5*dot(x,A*x)-dot(b,x);     % energie

x_CG=zeros(n,iters+1);              % aproximace reseni
r=b-A*x_CG(:,1);                    % reziduum
p=r;                                % pocatecni smer

fprintf('k = %2d, ||r|| = %.4e, E = %.12e\n', ...
    0,norm(r),E(x_CG(:,1)));

for k=1:iters
    Ap=A*p;
    alpha=dot(r,r)/dot(p,Ap);       % delka kroku
    x_CG(:,k+1)=x_CG(:,k)+alpha*p;
    r_new=r-alpha*Ap;               % aktualizace rezidua

    fprintf('k = %2d, ||r|| = %.4e, E = %.12e\n', ...
        k,norm(r_new),E(x_CG(:,k+1)));

    if norm(r_new)<tol
        break
    end

    beta=dot(r_new,r_new)/dot(r,r); % koeficient sdruzeni
    p=r_new+beta*p;                 % novy sdruzeny smer
    r=r_new;
end