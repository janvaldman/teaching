%% Pruseciky tri rovin: vyber jedne geometricke situace
% Kazdy radek matice A urcuje jednu rovinu: A(i,:)*x = b(i).
%
% VOLBA PRIKLADU:
%   variant = 1 ... roviny nemaji spolecny bod (metoda nejmensich ctvercu)
%   variant = 2 ... roviny se protinaji v jednom bode
%   variant = 3 ... roviny maji spolecnou prusecnici (primku)

clear; close all; clc;

variant = 1;                 % ZMENIT NA 2 NEBO 3

% Data tri prikladu
A{1} = [1  1  1;             % Neslucitelna soustava
        1 -1  2;
        2  0  3];
b{1} = [3; 2; 1];

A{2} = [1  1  1;             % Jednoznacne reseni
        1 -1  2;
        0  1  1];
b{2} = [3; 2; 2];

A{3} = A{1};                 % Nekonecne mnoho reseni
b{3} = [3; 2; 5];

% Kontrola volby
if ~ismember(variant,1:3)
    error('Promenna variant musi byt 1, 2 nebo 3.');
end

% Vybereme pouze jednu variantu
M = A{variant};
rhs = b{variant};

figure('Name',sprintf('Pruseciky tri rovin - varianta %d',variant));
hold on;
plotPlanes(M,rhs);

% Jednotny vypocet reseni pomoci Mooreovy-Penroseovy pseudoinverze
x = pinv(M)*rhs;

% Vzdy vykreslime prave jeden bod
plot3(x(1),x(2),x(3),'rp','MarkerSize',13,...
      'MarkerFaceColor','r','DisplayName','Reseni pinv');

title(sprintf('Reseni pomoci pseudoinverze: (%.3g, %.3g, %.3g)',x));
legend('Location','bestoutside');

%% Vykresleni pruseciku rovin s krychli [-10,10]^3
function plotPlanes(A,b)
    lim = [-10 10];
    colors = lines(3);

    % Vrcholy krychle
    [X,Y,Z] = ndgrid(lim,lim,lim);
    V = [X(:),Y(:),Z(:)];

    % Dvanact hran krychle (indexy vrcholu v matici V)
    E = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; ...
         3 7; 4 8; 5 6; 5 7; 6 8; 7 8];

    for i = 1:3
        a = A(i,:);
        c = b(i);

        % Pruseciky roviny se vsemi hranami krychle
        P = [];
        for e = 1:size(E,1)
            p1 = V(E(e,1),:);
            p2 = V(E(e,2),:);
            f1 = dot(a,p1)-c;
            f2 = dot(a,p2)-c;

            if abs(f1) < 1e-12
                P = [P; p1]; %#ok<AGROW>
            end
            if f1*f2 < 0
                t = f1/(f1-f2);
                P = [P; p1+t*(p2-p1)]; %#ok<AGROW>
            end
        end

        % Odstraneni duplicitnich bodu
        if isempty(P)
            continue
        end
        P = unique(round(P,12),'rows');

        % Serazeni vrcholu polygonu v rovine
        if size(P,1) >= 3
            center = mean(P,1);
            n = a(:)/norm(a);
            [~,idx] = min(abs(n));
            ref = zeros(3,1);
            ref(idx) = 1;
            u = cross(n,ref);
            u = u/norm(u);
            v = cross(n,u);

            D = P-center;
            angle = atan2(D*v,D*u);
            [~,order] = sort(angle);
            P = P(order,:);

            patch('Vertices',P,'Faces',1:size(P,1),...
                  'FaceColor',colors(i,:),...
                  'FaceAlpha',0.35,...
                  'EdgeColor','none',...
                  'DisplayName',sprintf('Rovina %d',i));
        elseif size(P,1) == 2
            % Degenerovany prusecik: usecka na hrane krychle
            plot3(P(:,1),P(:,2),P(:,3),'-',...
                  'Color',colors(i,:),'LineWidth',2,...
                  'DisplayName',sprintf('Rovina %d',i));
        end
    end

    xlabel('x'); ylabel('y'); zlabel('z');
    xlim(lim); ylim(lim); zlim(lim);
    axis equal;
    grid on;
    box on;
    view(3);
    set(gca,'Projection','perspective');
end
