n=15;                                             %maximalni pocet clenu
p=[];
for k=1:n
    p=[(-1)^(k-1)*(1/(2*k-1)), 0, p];             %koeficienty polynomu
    pi_k(k)=16*polyval(p,1/4)-4*polyval(p,79/401); %aproximace pi
end
e_k=pi_k-pi;                                      %absolutni chyba 
format long; 
fprintf('%2.0f:  %1.16f  %2.2d\n', [(1:n)' pi_k' abs(e_k)'].')
semilogy(abs(e_k),'o-');                          %graf v log skale

