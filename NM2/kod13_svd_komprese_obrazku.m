%% Komprese obrazku pomoci singularniho rozkladu (SVD)
clear; close all; clc;

% =========================================================
% NACTENI A PREVOD OBRAZKU
% =========================================================
inImage = imread('cb_market.jpeg');

% Prevod barevneho obrazu na odstiny sedi
if size(inImage,3) == 3
    inImage = rgb2gray(inImage);
end

A = double(inImage);

[m,n] = size(A);
originalElements = m*n;

fprintf('Rozmer obrazku: %d x %d\n',m,n);
fprintf('Pocet pixelu: %d\n',originalElements);

% =========================================================
% SINGULARNI ROZKLAD
% =========================================================
fprintf('\nPocitam singularni rozklad...\n');

[U,S,V] = svd(A,'econ');

sigma = diag(S);
r = length(sigma);

fprintf('Singularni rozklad dokoncen.\n');

% =========================================================
% TESTOVANE POCTY SINGULARNICH CISEL
% =========================================================
all_N = [5 20 80 320];
all_N = all_N(all_N <= r);

numN = length(all_N);

% =========================================================
% PRIPRAVA VYSLEDKU
% =========================================================
relativeError = zeros(numN,1);
storedElements = zeros(numN,1);
compressionRatio = zeros(numN,1);
memorySaving = zeros(numN,1);

% Celkova Frobeniova norma matice
normA = norm(sigma);

% =========================================================
% GRAFICKE OKNO
% =========================================================
figure('Color','w','Position',[100 100 1500 850]);

tl = tiledlayout(2,3,...
    'TileSpacing','compact',...
    'Padding','compact');

% =========================================================
% REKONSTRUKCE OBRAZKU
% =========================================================
for i = 1:numN

    N = all_N(i);

    fprintf('Rekonstruuji obrazek pro N = %d...\n',N);

    % Redukovany singularni rozklad
    Ured = U(:,1:N);
    Sred = S(1:N,1:N);
    Vred = V(:,1:N);

    % Rekonstrukce obrazku
    D = Ured*Sred*Vred';

    % Relativni chyba vypoctena ze singularnich cisel
    relativeError(i) = ...
        sqrt(sum(sigma(N+1:end).^2))/normA;

    % Pocet ulozenych prvku redukovaneho rozkladu
    storedElements(i) = N*(m+n+1);

    % Kompresni pomer a uspora pameti
    compressionRatio(i) = originalElements/storedElements(i);

    memorySaving(i) = ...
        100*(1-storedElements(i)/originalElements);

    % Zobrazeni rekonstruovaneho obrazku
    nexttile(i);
    imshow(uint8(D));

    title(sprintf('%d singularnich cisel',N));

end

% =========================================================
% GRAF RELATIVNI CHYBY
% =========================================================
ax1 = nexttile(5);
hold(ax1,'on');

semilogx(ax1,all_N,relativeError,'-o',...
    'LineWidth',1.8,...
    'MarkerSize',7);

% Popisky hodnot relativni chyby
for i = 1:numN

    text(ax1,all_N(i)*1.08,relativeError(i),...
        sprintf('%.3f',relativeError(i)),...
        'FontSize',10,...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','left');

end

grid(ax1,'on');
axis(ax1,'tight');

xlabel(ax1,'Pocet singularnich cisel');
ylabel(ax1,'Relativni chyba');
title(ax1,'Chyba aproximace');

% =========================================================
% GRAF POCTU ULOZENYCH PRVKU
% =========================================================
ax2 = nexttile(6);
hold(ax2,'on');

plot(ax2,all_N,storedElements/1e6,'-o',...
    'LineWidth',1.8,...
    'MarkerSize',7);

% Pocet prvku puvodniho obrazku
originalMillions = originalElements/1e6;

yline(ax2,originalMillions,'--r',...
    'LineWidth',1.5);

% Popisek puvodniho obrazku
text(ax2,8,originalMillions+0.04,...
    'Puvodni obraz',...
    'Color','r',...
    'FontSize',10,...
    'FontWeight','bold',...
    'VerticalAlignment','bottom');

% Nastaveni rozsahu os
xlim(ax2,[0 340]);
ylim(ax2,[0 originalMillions*1.08]);

grid(ax2,'on');

xlabel(ax2,'Pocet singularnich cisel');
ylabel(ax2,'Pocet ulozenych prvku (miliony)');
title(ax2,'Velikost redukovaneho rozkladu');

% =========================================================
% VYPIS VYSLEDKU
% =========================================================
fprintf('\nVYSLEDKY KOMPRESE OBRAZKU\n');
fprintf('================================================================\n');
fprintf('   N     Rel. chyba   Ulozene prvky   Komprese    Uspora pameti\n');
fprintf('----------------------------------------------------------------\n');

for i = 1:numN

    fprintf('%4d     %8.5f      %10d      %6.2f : 1      %7.2f %%\n',...
        all_N(i),...
        relativeError(i),...
        storedElements(i),...
        compressionRatio(i),...
        memorySaving(i));

end

fprintf('================================================================\n');
fprintf('Pocet prvku puvodniho obrazku: %d\n',originalElements);

% =========================================================
% EXPORT GRAFU
% =========================================================
fprintf('\nUkladam obrazek...\n');

drawnow;
print(gcf,'svd_komprese_obrazku.png','-dpng','-r300');

fprintf('Obrazek ulozen.\n');
fprintf('Program dokoncen.\n');