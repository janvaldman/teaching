close all

rozmer = 2;                 % rozměr obrázku: 2 nebo 3

bodX = [1 1 1]';            % bod mimo přímku
bodA = [-1/2 0 0]';         % první bod přímky
bodB = [1 2 2]';            % druhý bod přímky

% Při volbě 2D odstraníme třetí souřadnici všech bodů
if rozmer == 2
    bodX(end) = [];
    bodA(end) = [];
    bodB(end) = [];
end

vektorAB = bodB - bodA;     % směrový vektor přímky
vektorAX = bodX - bodA;     % vektor z bodu A do bodu X

% Výpočet projekce pomocí skalárního parametru
parametr = dot(vektorAX, vektorAB) / dot(vektorAB, vektorAB);
bodProjekce = bodA + parametr * vektorAB;

% Výpočet projekce pomocí projekční matice
maticeProjekce = (vektorAB * vektorAB') / dot(vektorAB, vektorAB);
bodProjekce2 = maticeProjekce * vektorAX + bodA;

hold on
oznac_bod(bodA, 'A');
oznac_bod(bodB, 'B');
oznac_bod(bodX, 'X');
oznac_bod(bodProjekce, 'P');

vykresli_usecku(bodA, bodB, '-go');        % přímka určená body A a B
vykresli_usecku(bodA, bodX, '-rx');        % spojnice bodů A a X
vykresli_usecku(bodProjekce, bodX, '--r'); % kolmá složka

axis image
grid on

if rozmer == 3
    view(3)                 % nastavení trojrozměrného pohledu
end

% Ověření kolmosti: výsledky by měly být přibližně nulové
format long
kontrolaKolmosti1 = dot(bodX - bodProjekce, vektorAB)
kontrolaKolmosti2 = dot(bodX - bodProjekce2, vektorAB)

function oznac_bod(bod, popisek)
if numel(bod) == 2
    text(bod(1), bod(2), popisek, 'FontSize', 20)
elseif numel(bod) == 3
    text(bod(1), bod(2), bod(3), popisek, 'FontSize', 20)
end
end

function vykresli_usecku(bod1, bod2, styl)
if numel(bod1) == 2
    plot([bod1(1) bod2(1)], [bod1(2) bod2(2)], styl)
elseif numel(bod1) == 3
    plot3([bod1(1) bod2(1)], ...
        [bod1(2) bod2(2)], ...
        [bod1(3) bod2(3)], styl)
end
end