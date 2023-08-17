%example 2.2, page 
A=[1 1 1; 1 -1 2; 2 0 3]; b=[3; 2; 1];      % no solution
%A=[1 1 1; 1 -1 2; 0 1 1]; b=[3; 2; 2];     % one solution
%A=[1 1 1; 1 -1 2; 2 0 3]; b=[3; 2; 5];      % infinitely many solutions

% overdetermined system - three lines intersection (wiki)
%A=[-2 -1; 3 -1; 1 -1]; b=[1; 2; -1];        % no solution
%A=[-2 -1; 1 -1; 2 -1]; b=[-1; -1; -1];      % one solution
%A=[1 -1; 1 -1; 1 -1]; b=[0; 0; 0];          % infinitely many solutions

format rat
%format long

[k,l]=size(A);   % dimensions of a matrix

Ab=[A b];      % auxiliary matrix
Ab
    
[AbEch,p]=rref(Ab);  % reduced row echelon form
AbEch

%  solution type
if ismember(l+1,p)
    fprintf('No solution! \n')
elseif nnz(ismember(1:l,p))==l
    fprintf('One solution! \n')
else
    fprintf('Infinitely many solutions! \n') 
end    

if 1
[AEch,p]=rref(A)            %reduced echelon form, p - pivots
AEchExt=-eye(l);    %minus 1 trick
AEchExt(p,:)=AEch(p,:)

n=AEchExt; n(:,p)=[]        %null space of A
n2=null(A);                 %another null space using Matlab 
%n./n2

check=norm(A*n)                         %check of null space

end