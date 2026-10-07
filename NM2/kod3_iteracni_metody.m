n=4; iters=15;                         % velikost matice a pocet iteraci 
A =(n+1)^2*gallery('tridiag',n,-1,2,-1);        % sestaveni ridke matice 
b=ones(n,1);
D=diag(diag(A)); L=tril(A,-1); U=triu(A,1);     % rozklad A=L+D+U
C_J=-D\(L+U); d_J=D\b;                          % Jacobi
C_GS=-(L+D)\U; d_GS=(L+D)\b;                    % Gauss-Seidel 
x_J=zeros(n,iters+1);  x_GS=x_J;                % iterace 
for k=1:iters
    x_J(:,k+1)=C_J*x_J(:,k)+d_J;  
    x_GS(:,k+1)=C_GS*x_GS(:,k)+d_GS; 
end
rho_J=max(abs(eigs(C_J)))                      % spektralni polomer 
rho_GS=max(abs(eigs(C_GS)))                    % spektralni polomer 
figure; 
subplot(1,2,1); plot(x_J,'b--'); axis square; 
subplot(1,2,2); plot(x_GS,'r:'); axis square; 

% Formatovani grafu
ymin = min([x_J(:); x_GS(:)]);
ymax = max([x_J(:); x_GS(:)]);
dy = 0.05*(ymax-ymin);

subplot(1,2,1);
xlim([0.75 n+0.25]); ylim([ymin-dy ymax+dy]);

subplot(1,2,2);
xlim([0.75 n+0.25]); ylim([ymin-dy ymax+dy]);