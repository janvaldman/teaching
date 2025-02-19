format rational
close all
x=[6; 0; 0];
%x=[1; 1; 1.2]
B=[1 0; 1 1; 1 2];    %vector generation the projection space
b1=B(:,1); b2=B(:,2);

%one way
lambda=(B'*B)\(B'*x)
x_b=B*lambda  % projection 

%another  
P=B*inv(B'*B)*B'   %projection matrix
x_b2=P*x

hold on
plot3([0 x(1)],[0 x(2)],[0 x(3)],'-rx')  % x vector 
plot3([0 b1(1)],[0,b1(2)],[0,b1(3)],'-go')   % b1 vector 
plot3([0 b2(1)],[0,b2(2)],[0,b2(3)],'-go')   % b2 vector 
plot3([0 x_b(1)],[0,x_b(2)],[0,x_b(3)],'-ro') % x_b vector 
plot3([x_b(1),x(1)],[x_b(2),x(2)],[x_b(3),x(3)],'--r') % orthogonal complement
%axis image
axis equal
%grid on
%enlarge_axis(0.1,0.1)
view(3)

%check of orthogonality
format long
check=B'*(x-x_b)
check2=B'*(x-x_b2)
