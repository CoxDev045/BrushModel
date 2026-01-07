clear;clc
syms h w z real

lambda = -z * w + 1i * 2 * w * sqrt(1 - z^2);
expanded_expr = 1 + h * lambda + 1/2 * (h * lambda)^2 + 1/6 * (h * lambda)^3 + 1/24 * (h * lambda)^4;
expanded_expr = simplify(expanded_expr)

expr = abs(expanded_expr) == 1;
solve(expr, h, ReturnConditions=true)
%%
% Stability region comparison for RK4 vs Taylor-5

% Complex plane grid
re = linspace(-5, 2, 1000);
im = linspace(-4, 4, 1000);
[Re, Im] = meshgrid(re, im);
Z = Re + 1i*Im;

R_ExplicitEuler = 1 + Z;
R_T2 = R_ExplicitEuler + Z.^2/2;
R_T3 = R_T2 + Z.^3/4;
% RK4 stability function
R_RK4_38 = R_T3 + Z.^4/24;
R_RK4 = R_T2 + Z.^3/6 + Z.^4/24;

% % 5th-order Taylor polynomial of exp(z)
% R_T5 = R_RK4 + Z.^5/120;

% Stability conditions
stab_ExplicitEuler = abs(R_ExplicitEuler) < 1;
stab_T2 = abs(R_T2) < 1;
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
