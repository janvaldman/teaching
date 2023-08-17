n=4;                                                        %pocet skupin
k=10;                                                       %pocet minci   
items=nchoosek(1:n+k-1,n-1);
numb=size(items,1); 
dist=diff([zeros(numb,1), items, (n+k)*ones(numb,1)],1,2)-1