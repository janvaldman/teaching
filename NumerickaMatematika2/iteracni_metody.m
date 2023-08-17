iters=15; n=4;                         %pocet iteraci a velikost matice
A=(n+1)^2*(2*diag(ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1));
b=ones(n,1);
D=diag(diag(A)); L=tril(A,-1); U=triu(A,1);     %rozklad A=L+D+U
C_Jac=-D\(L+U); d_Jac=D\b;                      %Jacobi
C_GS=-(L+D)\U; d_GS=(L+D)\b;                    %Gauss-Seidel 
x_Jac=zeros(n,iters+1);  x_GS=x_Jac;            %matice iteraci 
for k=1:iters
    x_Jac(:,k+1)=C_Jac*x_Jac(:,k)+d_Jac;      
    x_GS(:,k+1)=C_GS*x_GS(:,k)+d_GS;    
end
figure; 
subplot(1,2,1); plot(x_Jac,'b--'); axis square; 
subplot(1,2,2); plot(x_GS,'r:'); axis square; 
rho_Jac=max(eig(C_Jac))                    %spektralni polomer 
rho_GS=max(eig(C_GS))                      %spektralni polomer 