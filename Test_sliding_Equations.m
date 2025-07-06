clearvars; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%
ax(1) = 0;
ay(1) = 0;

vx(1) = 0;
vy(1) = 0;

dx(1) = 0;
dy(1) = 0;

tx(1) = 0;
ty(1) = 0;

vs(1) = 1e-16;
% omega = 0;
re = 0.2;
mu = 1;

time_initial = 0;
time_final = 0.2;
fs = 1000;
dt = 1 / fs;
N = time_final * fs;

time = linspace(time_initial, time_final, N);
p = 1.0;
omega = (linspace(0, 10, N));
omega_z = zeros(size(omega)); %linspace(1e-16, 0.1, N);


k = 0.27;
c = 1.4e-4;
m = 7.64e-10;

vrx = @(omega, omega_z, y, dy) omega * re + omega_z * (y + dy);
vry = @(omega, omega_z, x, dx) - omega_z * (x + dx);

for i = 1:N-1
    %%%%%%%%%%%% Calculate angles from initial conditions %%%%%%%%%%%%%%%%%%%%

    theta1(i) = atan( ( vry(omega(i), omega_z(i), 0, dx(i)) + vy(i) ) / (vx(i) - vrx(omega(i), omega_z(i), 0, dy(i)) ) );
    theta2(i) = theta1(i) - pi;
    vs(i) = ( vx(i) - vrx(omega(i), omega_z(i), 0, dy(i)) ) / cos(theta1(i));

    %%%%%%%%%%%%%%%%%%%% Calculate next iterations stresses %%%%%%%%%%%%%%%%%%
    
    tx(i+1) = mu * p * cos(theta2(i));
    ty(i+1) = mu * p * sin(theta2(i));
    
    %%%%%%%%%%%%%%%%% Solve ODE for n = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ax(i+1) = 1/m * (tx(i+1) - k * dx(i) - c * vx(i));
    ay(i+1) = 1/m * (ty(i+1) - k * dy(i) - c * vy(i));
    
    %%%% Calculate velocity for n = 1 %%%%%%%%%%%%
    
    vx(i+1) = vx(i) + ax(i+1) * dt;
    vy(i+1) = vy(i) + ay(i+1) * dt;
    
    %%%% Calculate displacement for n = 1 %%%%%%%%%%%%

    dx(i+1) = dx(i) + vx(i+1) * dt;
    dy(i+1) = dy(i) + vy(i+1) * dt;

end

close all;

figure
plot(time(2:N), vs)
hold on
plot(time, omega)
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend("Sliding velocity", "Angular Velocity")

figure
subplot(131)
plot(dx, dy)
grid on
xlabel('x')
ylabel('y')
title("Deformation")

subplot(132)
plot(vx, vy)
grid on
xlabel('x')
ylabel('y')
title("Deformation velocity")

subplot(133)
plot(ax, ay)
grid on
xlabel('x')
ylabel('y')
title("Deformation acceleration")

%%
clearvars; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')

time_initial = 0;
time_final = 0.25;
fs = 1e7;
dt = 1 / fs;
N = time_final * fs;

time = linspace(time_initial, time_final, N);
p = 0.274539158098840;% * (sin(2 * pi * time) + 1.1);
omega = (linspace(0.730730730730731, 10, N));
omega_z = zeros(size(omega)); %linspace(1e-16, 0.1, N);


for i = 1:N
    X_vec(:, i) = solve_brush_ode(omega(i), omega_z(i), 0, 0, p, dt);
end
%
figure
plot(time, X_vec(5, :))
hold on
plot(time, omega)
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend("Sliding velocity", "Angular Velocity")

figure
subplot(131)
plot(X_vec(6, :), X_vec(7, :))
grid on
xlabel('x')
ylabel('y')
title("Deformation")

subplot(132)
plot(X_vec(3, :), X_vec(4, :))
grid on
xlabel('x')
ylabel('y')
title("Deformation velocity")

subplot(133)
plot(X_vec(1, :), X_vec(2, :))
grid on
xlabel('x')
ylabel('y')
title("Deformation acceleration")

%%
function [X_vec] = solve_brush_ode(omega, omega_z, v0, alpha, p, dt)
    arguments
        omega   (1, 1) double
        omega_z (1, 1) double
        v0      (1, 1) double
        alpha   (1, 1) double
        p       (1, 1) double
        dt      (1, 1) double
    end

    % Inputs:
    % vrx, vry - initial velocities in x and y directions
    % mu, p - friction coefficient and pressure
    % m - mass
    % k_x, k_y - stiffness constants
    % c_x, c_y - damping coefficients
    % x, y - positions
    persistent ax ay vx vy delta_x delta_y vrx vry vs theta_1 theta_2 tauX tauY ;
    persistent m re kx ky cx cy mu y x;

    if isempty(ax)

        m = 7.64e-10; re = 30;
        kx = 0.689725460841901; ky = 0.369725460841901;
        cx = 1.78e-07; cy = 1.4e-04;
        mu = 0.022065759598109;
        x = -6.888160351067128; y = -7;

        ax = 0;
        ay = 0;
        
        vx = 0;
        vy = 0;

        delta_x = 0.008030803080308;
        delta_y = 0;

        vrx = 0.217048731900217;
        vry = 0;

        theta_1 = 0;
        theta_2 = -pi;

        vs = 0;

        tauX = 0.005539090261749;
        tauY = 0;
    end
    
    %%%%%%%%%%%% Calculate angles from initial conditions %%%%%%%%%%%%%%%%%%%%
    vrx = omega * re + omega_z * (y + delta_y) + v0 * cos(alpha);
    vry = - omega_z * (x + delta_x) + v0 * sin(alpha);
    
    theta_1 = atan( ( vy - vry ) / (vx - vrx) );
    theta_2 = theta_1 - pi;
    vs = ( vx - vrx ) / cos(theta_1);

    %%%%%%%%%%%%%%%%%%%% Calculate next iterations stresses %%%%%%%%%%%%%%%%%%

    tauX = mu * p * cos(theta_2);
    tauY = mu * p * sin(theta_2);
    
    % % % %%%%%%%%%%%%%%%%% Solve ODE for n = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % 
    % % % ax = 1/m * (tauX - kx * delta_x - cx * vx);
    % % % ay = 1/m * (tauY - ky * delta_y - cy * vy);
    % % % 
    % % % %%%% Calculate velocity for n = 1 %%%%%%%%%%%%
    % % % 
    % % % vx = vx + ax * dt;
    % % % vy = vy + ay * dt;
    % % % 
    % % % %%%% Calculate displacement for n = 1 %%%%%%%%%%%%
    % % % 
    % % % delta_x = delta_x + vx * dt;
    % % % delta_x = delta_x + vy * dt;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%% ChatGPT Runge-Kutta %%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % % Define system of ODEs for RK4
    % function dX = brush_dynamics(~, X)
    %     dx = X(3);
    %     dy = X(4);
    %     dvx = (1/m) * (tauX - kx * X(1) - cx * X(3));
    %     dvy = (1/m) * (tauY - ky * X(2) - cy * X(4));
    %     dX = [dx; dy; dvx; dvy];
    % end
    % 
    % % RK4 Integration
    % X = [delta_x; delta_y; vx; vy];  % Current state
    % k1 = dt * brush_dynamics(0, X);
    % k2 = dt * brush_dynamics(0, X + 0.5 * k1);
    % k3 = dt * brush_dynamics(0, X + 0.5 * k2);
    % k4 = dt * brush_dynamics(0, X + k3);
    % 
    % X_next = X + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    % 
    % % Update state variables
    % delta_x = X_next(1);
    % delta_y = X_next(2);
    % vx = X_next(3);
    % vy = X_next(4);

    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%%%%%%%%%%%%%%%%%%%%%% Verlet Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 
    % % vx = vx + 0.5 * ax * dt;
    % % vy = vy + 0.5 * ay * dt;
    % % 
    % % delta_x = delta_x + vx * dt + 0.5 * ax * dt^2;
    % % delta_y = delta_y + vy * dt + 0.5 * ay * dt^2;
    % % 
    % % ax_new = 1/m * (tauX - kx * delta_x - cx * vx);
    % % ay_new = 1/m * (tauY - ky * delta_y - cy * vy);
    % % 
    % % vx = vx + 0.5 * (ax + ax_new) * dt;
    % % vy = vy + 0.5 * (ay + ay_new) * dt;

    % % Update values for n-1
    % prev1_delta_x = delta_x;
    % prev1_delta_y = delta_y;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% 4th Order Yoshida %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w0 = -(2^(1/3)) / (2 - (2^(1/3)));
    w1 = -(1) / (2 - (2^(1/3)));
    
    c1 = w1 / 2; c4 = c1;
    c2 = (w0 + w1) / 2 ;c3 = c2;

    d1 = w1; d3 = d1;
    d2 = w0;

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta_x = delta_x + c1 .* vx .* dt;
    delta_y = delta_y + c1 .* vy .* dt;

    % Update acceleration with temporary values
    ax = (tauX - kx .* delta_x - cx .* vx) ./ m;
    ay = (tauY - ky .* delta_y - cy .* vy) ./ m;
    
    % Update temporary velocity with new accel
    vx = vx + d1 .* ax .* dt;
    vy = vy + d1 .* ay .* dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta_x = delta_x + c2 .* vx .* dt;
    delta_y = delta_y + c2 .* vy .* dt;

    % Update acceleration with temporary values
    ax = (tauX - kx .* delta_x - cx .* vx) ./ m;
    ay = (tauY - ky .* delta_y - cy .* vy) ./ m;
    
    % Update temporary velocity with new accel
    vx = vx + d2 .* ax .* dt;
    vy = vy + d2 .* ay .* dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta_x = delta_x + c3 .* vx .* dt;
    delta_y = delta_y + c3 .* vy .* dt;

    % Update acceleration with temporary values
    ax = (tauX - kx .* delta_x - cx .* vx) ./ m;
    ay = (tauY - ky .* delta_y - cy .* vy) ./ m;
    
    % Update temporary velocity with new accel
    vx = vx + d3 .* ax .* dt;
    vy = vy + d3 .* ay .* dt;

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update brush displacement with correctors
    delta_x = delta_x + c4 .* vx .* dt;
    delta_y = delta_y + c4 .* vy .* dt;

    X_vec = [ax; %1
         ay;     %2
         vx;     %3
         vy;     %4
         vs;     %5
         delta_x;%6
         delta_y;%7
         vrx;    %8
         vry;    %9
         theta_1;%10
         theta_2;%11
         tauX;   %12
         tauY];  %13
end