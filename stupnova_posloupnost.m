%seq=[1, 2, 3, 4, 4];                         % jedna vstupni posloupnost
seq=[6, 6, 5, 4, 3, 3, 3, 3, 2, 1, 0];       % druha vstupni posloupnost
while ~isempty(seq)
    seq=sort(seq,'descend' );                % serazeni 
    fprintf('%d  ',seq); fprintf('\n',seq)   % vypis
    p=seq(1); 
    if (p<0) || (numel(seq)<p)  
        fprintf('Posloupnost neni grafova. \n'); return;    
    end
    if (p==0) && (numel(seq)==1)  
        fprintf('Posloupnost je grafova. \n'); return;  
    end
    seq(1)=[];                               % vymazani prvniho prvku
    seq(1:p)=seq(1:p)-1;                     % odecteni 1 od dalsi clenu
end