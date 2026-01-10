clear;clc
syms h w zeta real
% z = h * w;
% A = [(1 - z^2/2), (z/w - zeta * z^2/(2*w));
%     (-z*w +z^3*w/4 + z^3 * zeta / (4 * w) ) ,( 1 - z*zeta + (z/w)^2 - z^4 / (4*w) + z^2*zeta / (4*w) )];
% 
% lambda = eig(A)

A = [0, 1;
    -w^2, -zeta*w];

[V, D] = eig(A);
simplify(D)

% % lambda = -z * w + 1i * 2 * w * sqrt(1 - z^2);
% expanded_expr = 1 + h * lambda + 1/2 * (h * lambda)^2 + 1/6 * (h * lambda)^3 + 1/24 * (h * lambda)^4;
% expanded_expr = simplify(expanded_expr)
% 
% expr = abs(expanded_expr) == 1;
% solve(expr, h, ReturnConditions=true)
%%
% Stability region comparison for RK4 vs Taylor-5

% Complex plane grid
re = linspace(-5, 2, 1000);
im = linspace(-4, 4, 1000);
[Re, Im] = meshgrid(re, im);
Z = Re + 1i*Im;

% stability function
R_ExplicitEuler = 1 + Z;
R_ImplicitEuler = 1./(1 + Z);
R_ImplicitTrapz = (1 + Z/2) ./ (1 - Z/2);
R_RK4 = 1 + Z + Z.^2/2 + Z.^3/6 + Z.^4/24;

% Stability conditions
stab_ExplicitEuler = abs(R_ExplicitEuler) < 1;
stab_T2 = abs(R_ImplicitTrapz) < 1;
stab_T3 = abs(R_T3) < 1;

stab_RK4 = abs(R_RK4) < 1;
stab_RK4_38 = abs(R_RK4_38) < 1;
% stab_T5  = abs(R_T5)  < 1;

% Plot
figure;
hold on;
contour(Re, Im, stab_ExplicitEuler, [1 1], 'LineWidth', 2); 
contour(Re, Im, stab_T2, [1 1], 'LineWidth', 2); 
contour(Re, Im, stab_T3, [1 1], 'LineWidth', 2); 

contour(Re, Im, stab_RK4, [1 1], 'LineWidth', 2); 
contour(Re, Im, stab_RK4_38,  [1 1], '--', 'LineWidth', 2);
hold off
axis equal;
xlabel('Re(z)');
ylabel('Im(z)');
legend('ExplicitEuler', 'Taylor-2', 'Taylor-3', 'RK4', 'Taylor-5');
title('Stability Regions');
grid on;
%%

% Stability region for leapfrog with damping

omega = 1;                     % normalize frequency
h_vals = linspace(0, 3, 1000);  % h*omega range
zeta_vals = linspace(0, 2.1, 1000);

[H, Z] = meshgrid(h_vals, zeta_vals);
stable = false(size(H));
% RK4 stability function
R = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;

for i = 1:size(H,1)
    for j = 1:size(H,2)

        h = H(i,j);
        zeta = Z(i,j);

        % % Leapfrog update matrix
        % a = 1 - (h^2*omega^2)/2;
        % b = h*(1 - h*zeta*omega/2);
        % c = (h^3*omega^4)/4 + (zeta*h^3*omega^2)/4 - h*omega^2;
        % d = h^2 - (h^4*omega^3)/4 + (h^2*omega*zeta)/4 - h*omega*zeta + 1;

        % M = [a, b;
        %      c, d];
        % 
        % eigvals = eig(M);

        % if max(abs(eigvals)) <= 1
        %     stable(i,j) = true;
        % end

        % Continuous-time eigenvalues
        disc = 1 - zeta^2/4;

        if disc >= 0
            lam1 = -zeta*omega/2 + 1i*omega*sqrt(disc);
            lam2 = conj(lam1);
        else
            lam1 = -zeta*omega/2 + omega*sqrt(-disc);
            lam2 = -zeta*omega/2 - omega*sqrt(-disc);
        end

        if max(abs([R(h*lam1)* R(h*lam2)])) <= 1
            stable(i,j) = true;
        end
    end
end

% Plot stability region
figure;
imagesc(h_vals, zeta_vals, stable);
set(gca,'YDir','normal');
xlabel('h \omega');
ylabel('\zeta');
title('Leapfrog Stability Region (Damped Oscillator)');
colorbar;
