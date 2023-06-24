n=32;                                           % pocet radek trojuhelnika
p=1;
A=zeros(n);
for i=1:n;
    A(i,1:i)=p; 
    p=mod([p 0]+[0 p],2);
end;
spy(A)                                          % vykresleni 