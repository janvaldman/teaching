close all
clear all

A=[1 0; 0 5];
A=[1 1; 1 5];
const=2;

[X,Y]=meshgrid(linspace(-const*1.1,const*1.1,100));
Z=X;
for i=1:size(X,1)
    for j=1:size(X,2)
        Z(i,j)=[X(i,j); Y(i,j)]'*A*[X(i,j); Y(i,j)];        
    end
end

figure
[c,h] = contour(X,Y,Z,'LevelStep',2,'Fill','on');
hold on
[c2,h2] = contour(X,Y,Z,[const const], 'Color','m');
set(h2,'LineStyle','-')
axis equal

[V, D]=eig(A)
v1=V(:,1); v2=V(:,2);
lambda1=D(1,1); lambda2=D(2,2);
ax1=sqrt(const/lambda1); ax2=sqrt(const/lambda2);
v1=ax1*v1; v2=ax2*v2;

plot([0 v1(1)]',[0 v1(2)]','y-');
plot([0 v2(1)]',[0 v2(2)]','y-');