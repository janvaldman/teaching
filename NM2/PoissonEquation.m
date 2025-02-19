close all
n = 20;                       
h = 1/(n + 1);
A = (1/h^2)*gallery('tridiag',n,-1,2,-1);        % sestaveni ridke matice 
x = h*linspace(1,n,n)';
b = x.^2;
u = A \ b;
scatter(x,u,'b','MarkerFaceColor','b')
uExact=@(y)(y-y.^4)/12;
hold on
xx=linspace(0,1)
plot(xx,uExact(xx),'r-')
hold off
e=max(abs(u-uExact(x)))


