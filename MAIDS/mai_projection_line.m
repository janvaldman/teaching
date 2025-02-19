format rat
close all

X=[1 1 1]';     % point outside the line
A=[-0.5 0 0]'; B=[1 2 2]';     % line points 
A=[0 0 0]'; B=[1 2 2]';     % line points 

% reduction to 2D 
X(end)=[]; A(end)=[]; B(end)=[]; 

AB=B-A;    %vector generating the projection space
AX=X-A;
% one way
lambda=dot(AX,AB)/dot(AB,AB)
PX=AB*lambda+A  % projection 

% another way
P=(1/dot(AB,AB))*(AB*AB');   %projection matrix
PX2=P*AX+A

hold on
text_own(A,'A');
text_own(B,'B');
text_own(X,'X');
text_own(PX,'PX');

plot_own(A,B,'-go'); 
plot_own(A,X,'-rx'); 
plot_own(PX,X,'--r'); % % orthogonal complement

axis image
grid on
%enlarge_axis(0.1,0.1)
%view(3)

%check of orthogonality
format long
check_orthogonality=dot(X-PX,AB)
check2_orthogonality=dot(X-PX2,AB)

function text_own(A,string)
    if numel(A)==2
        text(A(1),A(2),string,'FontSize',20)
    elseif numel(A)==3
        text(A(1),A(2),A(3),string,'FontSize',20)
    end
end

function plot_own(A,B,format)
    if numel(A)==2
        plot([A(1) B(1)],[A(2),B(2)],format);
    elseif numel(A)==3
        plot3([A(1) B(1)],[A(2),B(2)],[A(3),B(3)],format);
    end
end






