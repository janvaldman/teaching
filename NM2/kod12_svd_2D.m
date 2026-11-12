% Singularni rozklad – geometricka interpretace ve 2D
clear; close all; clc;

% =========================================================
% VSTUPNI MATICE
% =========================================================
A = [2 1; 0 1];

% Singularni rozklad
[U,S,V] = svd(A);
sigma = diag(S);

v1 = V(:,1); v2 = V(:,2);
u1 = U(:,1); u2 = U(:,2);

w1 = A*v1;
w2 = A*v2;

% Jednotkova kruznice a jeji obraz
t = linspace(0,2*pi,500);
X = [cos(t); sin(t)];
Y = A*X;

% Barvy
col1 = [0.90 0.15 0.05];
col2 = [0.05 0.25 0.95];

%% =========================================================
% GRAFICKE OKNO
% =========================================================
figure('Color','w','Position',[100 100 1250 570]);

tl = tiledlayout(1,2,...
    'TileSpacing','none',...
    'Padding','compact');

% =========================================================
% LEVY GRAF: VSTUPNI JEDNOTKOVA KRUZNICE
% =========================================================
ax1 = nexttile;
hold(ax1,'on');
grid(ax1,'on');
box(ax1,'on');

plot(ax1,X(1,:),X(2,:),'k','LineWidth',2);

quiver(ax1,0,0,v1(1),v1(2),0,...
    'Color',col1,'LineWidth',2.5,'MaxHeadSize',0.25);

quiver(ax1,0,0,v2(1),v2(2),0,...
    'Color',col2,'LineWidth',2.5,'MaxHeadSize',0.25);

axis(ax1,'equal');
xlim(ax1,[-1.35 1.35]);
ylim(ax1,[-1.35 1.35]);

xlabel(ax1,'$x_1$','Interpreter','latex','FontSize',17);
ylabel(ax1,'$x_2$','Interpreter','latex','FontSize',17);

title(ax1,'\textbf{Input: unit circle}',...
    'Interpreter','latex','FontSize',18);

set(ax1,'FontSize',14,'FontName','Arial','LineWidth',0.8);

addVectorLabel(ax1,v1,'$v_1$',col1,0.08);
addVectorLabel(ax1,v2,'$v_2$',col2,0.08);

% =========================================================
% PRAVY GRAF: VYSTUPNI ELIPSA
% =========================================================
ax2 = nexttile;
hold(ax2,'on');
grid(ax2,'on');
box(ax2,'on');

plot(ax2,Y(1,:),Y(2,:),'k','LineWidth',2);

quiver(ax2,0,0,w1(1),w1(2),0,...
    'Color',col1,'LineWidth',2.5,'MaxHeadSize',0.18);

quiver(ax2,0,0,w2(1),w2(2),0,...
    'Color',col2,'LineWidth',2.5,'MaxHeadSize',0.18);

axis(ax2,'equal');

maxY = max(abs(Y(:)));
limY = 1.5*max(1,maxY);

xlim(ax2,[-limY limY]);
ylim(ax2,[-limY limY]);

xlabel(ax2,'$y_1$','Interpreter','latex','FontSize',17);
ylabel(ax2,'$y_2$','Interpreter','latex','FontSize',17);

title(ax2,'\textbf{Output: ellipse $Ax$}',...
    'Interpreter','latex','FontSize',18);

set(ax2,'FontSize',14,'FontName','Arial','LineWidth',0.8);

addVectorLabel(ax2,w1,...
    sprintf('$\\sigma_1 u_1\\;(\\sigma_1=%.2f)$',sigma(1)),...
    col1,0.12);

addVectorLabel(ax2,w2,...
    sprintf('$\\sigma_2 u_2\\;(\\sigma_2=%.2f)$',sigma(2)),...
    col2,0.12);

% =========================================================
% HLAVNI TITULEK
% =========================================================
sgtitle(tl,...
    '\textbf{Singular Value Decomposition: $A=U\Sigma V^T$}',...
    'Interpreter','latex','FontSize',23);

%% =========================================================
% NUMERICKE OVERENI
% =========================================================
fprintf('\nMatice A:\n');
disp(A);

fprintf('Singularni cisla:\n');
disp(sigma);

fprintf('Prave singularni vektory V:\n');
disp(V);

fprintf('Leve singularni vektory U:\n');
disp(U);

fprintf('Overeni vztahu A*V = U*S:\n');
disp(A*V-U*S);

fprintf('Chyba rekonstrukce:\n');
fprintf('||A-U*S*V''|| = %.3e\n',norm(A-U*S*V'));

%% =========================================================
% EXPORT GRAFU
% =========================================================
exportgraphics(gcf,'svd_geometry2D.png','Resolution',300);
exportgraphics(gcf,'svd_geometry2D.pdf','ContentType','vector');


%% =========================================================
% LOKALNI FUNKCE: UMISTENI POPISKU UVNITR GRAFU
% =========================================================
function addVectorLabel(ax,w,label,col,offset)

    if norm(w)<1e-12
        return;
    end

    d = w/norm(w);

    xl = xlim(ax);
    yl = ylim(ax);

    scale = max(diff(xl),diff(yl));
    q = w + offset*scale*d;

    if d(1)>=0
        ha = 'right';
    else
        ha = 'left';
    end

    if d(2)>=0
        va = 'bottom';
    else
        va = 'top';
    end

    h = text(ax,q(1),q(2),label,...
        'Interpreter','latex',...
        'Color',col,...
        'FontSize',16,...
        'FontWeight','bold',...
        'HorizontalAlignment',ha,...
        'VerticalAlignment',va,...
        'Clipping','on');

    drawnow;

    % Skutecne rozmery textoveho popisku
    ext = get(h,'Extent');

    % Okraje uvnitr os
    mx = 0.035*diff(xl);
    my = 0.035*diff(yl);

    % Korekce vodorovne polohy
    if ext(1)<xl(1)+mx
        dx = xl(1)+mx-ext(1);
    elseif ext(1)+ext(3)>xl(2)-mx
        dx = xl(2)-mx-ext(1)-ext(3);
    else
        dx = 0;
    end

    % Korekce svisle polohy
    if ext(2)<yl(1)+my
        dy = yl(1)+my-ext(2);
    elseif ext(2)+ext(4)>yl(2)-my
        dy = yl(2)-my-ext(2)-ext(4);
    else
        dy = 0;
    end

    pos = get(h,'Position');
    pos(1) = pos(1)+dx;
    pos(2) = pos(2)+dy;
    set(h,'Position',pos);

end

%exportgraphics(gcf, 'svd_geometry2D.png', 'Resolution', 300);
%exportgraphics(gcf, 'svd_geometry2D.pdf', 'ContentType', 'vector');