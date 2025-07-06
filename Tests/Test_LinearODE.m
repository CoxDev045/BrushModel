clear; close all; clc;

k = 1;
m = 1;
c = 1;

params = [m, c, k];

cc = 2 * sqrt(k * m);
zeta = c / cc;

fprintf('With paramerters [m, c, k] = [%.2f, %.2f, %.2f], we have damping ratio Zeta = %.2f \n', [m, c, k], zeta)

A = [0, 1;
    -k/m, -c/m];

[V, D] = eig(A)

X0 = [1; 0];
Tspan = [0, 10];
time = linspace(Tspan(1), Tspan(2), 1000);

for i = 1:length(time)
    exp_sol(:, i) = expm(A * time(i)) * X0;
    eig_sol(:, i) = (V * diag(exp(diag(D * time(i)))) / V) * X0;
end

% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'NormControl', 'on', ...
%                 'OutputFcn', @odeplot, 'Stats','on');
% 
% sol = ode45(@(t, x) systemDynamics(t, x, A), time, X0, options)


figure
subplot(211)
plot(time, exp_sol.', '-o')
grid on

subplot(212)
plot(time, eig_sol, '-o')
grid on
%%
clc; close all

t_ind = find(output.time >= 90 & output.time <= 160);
time = output.time(t_ind);

edge0 = 92.5 - time(1);
edge1 = 115 - time(1);
edge2 = 132 - time(1);
edge3 = 157 - time(1);

time = time - time(1);
disp = output.disp_x(t_ind);

data_range = max(disp) - min(disp);
y = data_range * smoothstep(edge0, edge1, time) .* (1 - smoothstep(edge2, edge3, time)) + min(output.disp_x);

v0 = gradient(y) ./ gradient(time);

figure
subplot(211)
plot(time, y)
hold on
plot(time, disp)
grid on
legend('Smootthstep')
hold off

subplot(212)
plot(time, v0)
% hold on
% plot(time, gradient(disp) ./ gradient(time))
grid on
legend('Smootthstep')

save('Experimental_V0.mat', "v0")
%%

function Xdot = systemDynamics(t, x, A)
    
    Xdot = A * [x(1); x(2)];
end
