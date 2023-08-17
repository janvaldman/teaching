A=[5 4 1; 5 5 0; 0 0 5; 1 0 4]

% decomposing the image using singular value decomposition
[U,S,V]=svd(A);

all_N=[1,2,3];

dispEr = [];
numSVals = [];
nonzerosTotal = [];

%figure;
for i=1:numel(all_N)
    N=all_N(i);
    
    Ureduced=U(:,1:N);
    Sreduced=S(1:N,1:N);
    Vreduced=V(:,1:N);
    
    % Construct an Image using the selected singular values
    D=Ureduced*Sreduced*Vreduced'

    error=sum(sum((A-D).^2));

    % store vals for display
    nonzerosTotal=[nonzerosTotal numel(Ureduced)+numel(Sreduced)+numel(Vreduced)];
    dispEr = [dispEr; error];
    numSVals = [numSVals; N];
    
end

% dislay the error graph
figure; 
subplot(2,1,1)
title('Error in compression');
plot(numSVals, dispEr,'-x');
grid on
xlabel('Number of Singular Values used');
ylabel('Error');

subplot(2,1,2)
title('Nonzeros in matrices');
plot(numSVals, nonzerosTotal,'-x');
grid on
xlabel('Number of Singular Values used');
ylabel('Nonzeros');

