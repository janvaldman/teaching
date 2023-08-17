inImage=imread('cb_market.jpeg');
%reading and converting the image
inImage=rgb2gray(inImage);
inImageD=double(inImage);
figure;
buffer = sprintf('original picture');
imshow(uint8(inImageD));
title(buffer);

% decomposing the image using singular value decomposition
[U,S,V]=svd(inImageD);
Sdiag=diag(S);

% Using different number of singular values (diagonal of S) to compress and
% reconstruct the image
dispEr = [];
numSVals = [];
nonzerosTotal = [];
all_N=5:25:300;


all_N=ceil(linspace(5,min(size(S)),4));
all_N=[5 20 80 320];

%figure;
for i=1:numel(all_N)
    N=all_N(i);
    
    Ureduced=U(:,1:N);
    Sreduced=S(1:N,1:N);
    Vreduced=V(:,1:N);
    Sdiagreduced=Sdiag(1:N);
    
    % Construct an Image using the selected singular values
    %D=Ureduced*Sreduced*Vreduced';

    D=Ureduced*diag(Sdiagreduced)*Vreduced';

    % display and compute error
    figure;
    %subplot(ceil(sqrt(numel(all_N))),ceil(sqrt(numel(all_N))),i)
    buffer = sprintf('%d singular values', N);
    imshow(uint8(D));
    title(buffer);
    error=sum(sum((inImageD-D).^2));

    % store vals for display
    nonzerosTotal=[nonzerosTotal numel(Ureduced)+numel(Sdiagreduced)+numel(Vreduced)];
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

