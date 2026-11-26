%% Nelinearni kyvadlo jako soustava ODR prvniho radu
%   phi'' + (g/l)*sin(phi) = 0,   phi(0) = pi/2,  phi'(0) = 1,  t in [0, 20]
% Substituce omega = phi' vede na soustavu pro y = (phi, omega):
%   phi'   = omega
%   omega' = -(g/l)*sin(phi)
% Reseni RK4 s pevnym krokem porovname s vestavenym resicem ode45.
clear; close all; clc;

g = 9.81; l = 1;                              % tihove zrychleni, delka kyvadla
F = @(t,y) [y(2); -(g/l)*sin(y(1))];          % prava strana soustavy
y0 = [pi/2; 1]; t0 = 0; tMax = 20;

% RK4 s pevnym krokem
n = 2000;
[t, Y] = rk4_soustava(F, t0, tMax, y0, n);

% Referencni reseni ode45 s prisnou toleranci
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[tRef, YRef] = ode45(F, [t0 tMax], y0, opts);
phiRef = interp1(tRef, YRef(:,1), t, 'spline');
fprintf('RK4 (n = %d): max |phi - phi_ode45| = %.2e\n', n, max(abs(Y(:,1) - phiRef)));

% Celkova energie (na jednotku hmotnosti) by se mela zachovavat
E = @(Y) 0.5*l^2*Y(:,2).^2 - g*l*cos(Y(:,1));
fprintf('RK4: zmena energie max |E - E(0)| = %.2e\n', max(abs(E(Y) - E(Y(1,:)))));

figure;
subplot(2,1,1);
plot(t, Y(:,1), 'b-', tRef, YRef(:,1), 'r--');
xlabel('t'); ylabel('\phi(t)'); legend('RK4', 'ode45'); grid on;
title('Uhel kyvadla');
subplot(2,1,2);
plot(Y(:,1), Y(:,2), 'b-');
xlabel('\phi'); ylabel('\omega'); grid on; axis equal;
title('Fazovy portret');

%% RK4 pro soustavy (y je sloupcovy vektor, Y(i,:) je stav v case t(i))
function [t, Y] = rk4_soustava(F, t0, tMax, y0, n)
    t = linspace(t0, tMax, n+1)'; h = t(2) - t(1);
    Y = zeros(n+1, numel(y0)); Y(1,:) = y0(:)';
    for i = 1:n
        y  = Y(i,:)';
        k1 = F(t(i),       y);
        k2 = F(t(i) + h/2, y + h/2*k1);
        k3 = F(t(i) + h/2, y + h/2*k2);
        k4 = F(t(i) + h,   y + h*k3);
        Y(i+1,:) = (y + h/6*(k1 + 2*k2 + 2*k3 + k4))';
    end
end
