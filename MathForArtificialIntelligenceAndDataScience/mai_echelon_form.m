% example 2.2, page 20 
A=[1 1 1; 1 -1 2; 2 0 3]; b=[3; 2; 1];  % no solution
A=[1 1 1; 1 -1 2; 0 1 1]; b=[3; 2; 2];  % one solution
%A=[1 1 1; 1 -1 2; 2 0 3]; b=[3; 2; 5];  % infinitely many solutions

% example 2.6, page 29, infinitely many solutions 
A=[-2, 4, -2, -1, 4; 4, -8, 3, -3, 1; 1, -2, 1, -1, 1; 1, -2, 0, -3, 4];
a=-1; b=[-3; 2; 0; a]; 

% overdetermined system - three lines intersection (wiki)
%A=[-2 -1; 3 -1; 1 -1]; b=[1; 2; -1];        % no solution
%A=[-2 -1; 1 -1; 2 -1]; b=[-1; -1; -1];      % one solution
%A=[1 -1; 1 -1; 1 -1]; b=[1; 1; 1];          % infinitely many solutions
%A=[1 0 1 1; 0 1 0 2; 0 0 0 0]; b=[3; 1; 0]; % infinitely many solutions

format rat
%format long

[k,l]=size(A);   % dimensions of a matrix

% looking for the null space of A 
A
[AEch,jb]=rref(A)            %reduced echelon form, p - pivots
AEchExt=-eye(l);            %minus 1 trick
AEchExt(jb,:)=AEch(1:numel(jb),:)

n=AEchExt; n(:,jb)=[]        %null space of A
n2=null(A);                 %another null space using Matlab 
%n./n2
check=norm(A*n)             %check of null space

fprintf('-----\n')
% looking for the inverse of A 
if k==l
   AI=[A eye(k)];
   [AAI jb]=rref(AI);
   if numel(jb)==k
      Ainv=AAI(:,k+1:2*k) 
   end
end
fprintf('-----\n')

% looking for solvability
Ab=[A b]      % auxiliary matrix
[AbEch,jb]=rref(Ab)  % reduced row echelon form, p contains pivot collumns

%  solution type
if ismember(l+1,jb)
    fprintf('No solution! \n');
    type_sol=0;
elseif nnz(ismember(1:l,jb))==l
    fprintf('One solution! \n')
    type_sol=1;
else
    fprintf('Infinitely many solutions! \n') 
    type_sol=2;


end    