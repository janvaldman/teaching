% Nelinearni kyvadlo  phi'' + (g/l)*sin(phi) = 0,  phi(0) = pi/2, phi'(0) = 1
% Substituce omega = phi' vede na soustavu prvniho radu pro y = (phi, omega):
%   phi' = omega,   omega' = -(g/l)*sin(phi)
g = 9.81; l = 1;                        % tihove zrychleni, delka kyvadla
F = @(t,y) [y(2); -(g/l)*sin(y(1))];    % prava strana soustavy y' = F(t,y)
y0 = [pi/2; 1]; t0 = 0; tMax = 2;      % pocatecni podminka a konec intervalu
n = 2000; h = (tMax - t0)/n; t = t0 + h*(0:n);

Y = zeros(2,n+1); Y(:,1) = y0;          % RK4 pro soustavu: sloupec Y(:,i) je stav v case t(i)
for i = 1:n
    k1 = F(t(i),       Y(:,i));
    k2 = F(t(i) + h/2, Y(:,i) + h/2*k1);
    k3 = F(t(i) + h/2, Y(:,i) + h/2*k2);
    k4 = F(t(i+1),     Y(:,i) + h*k3);
    Y(:,i+1) = Y(:,i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

% Srovnani s vestavenym resicem ode45 v tychz casech t
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, Yref] = ode45(F, t, y0, opts);

% Energie na jednotku hmotnosti by se mela zachovavat
E = 0.5*l^2*Y(2,:).^2 - g*l*cos(Y(1,:));

fprintf('\n  %-33s = %.2e\n', 'krok h', h);
fprintf('  %-33s = %.2e\n', 'max |phi_RK4 - phi_ode45|', max(abs(Y(1,:) - Yref(:,1)')));
fprintf('  %-33s = %.2e\n\n', 'max |E - E(0)|', max(abs(E - E(1))));

% Perioda kmitu
% a) numericky: casy pruchodu phi nulou zdola nahoru (linearni interpolace)
phi = Y(1,:);
idx = find(phi(1:end-1) < 0 & phi(2:end) >= 0);
tNul = t(idx) - phi(idx).*(t(idx+1) - t(idx))./(phi(idx+1) - phi(idx));
T_num = mean(diff(tNul));
% b) presne: amplituda ze zakona zachovani energie a elipticky integral
phiMax = acos(cos(y0(1)) - l*y0(2)^2/(2*g));
T_ex = 4*sqrt(l/g)*ellipke(sin(phiMax/2)^2);
% c) linearizace pro male kmity (sin(phi) ~ phi): T_0 = 2*pi*sqrt(l/g)
T_0 = 2*pi*sqrt(l/g);

fprintf('  %-33s = %.4f rad\n', 'amplituda phi_max', phiMax);
fprintf('  %-33s = %.7f s\n', 'perioda z reseni RK4', T_num);
fprintf('  %-33s = %.7f s\n', 'perioda presna (ellipke)', T_ex);
fprintf('  %-33s = %.7f s\n\n', 'perioda malych kmitu T_0', T_0);

figure(1)                               % uhel v case
plot(t, Y(1,:), '-', t, Yref(:,1), '--', 'LineWidth', 1.5)
grid on; xlabel('t'); ylabel('\phi(t)')
legend('RK4', 'ode45', 'Location', 'southoutside', 'Orientation', 'horizontal')
title('Uhel kyvadla')

figure(2)                               % fazovy portret
plot(Y(1,:), Y(2,:), 'LineWidth', 1.5)
grid on; xlabel('\phi'); ylabel('\omega')
title('Fazovy portret')
