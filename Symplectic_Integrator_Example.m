clearvars; close all; clc;
set(0, "DefaultFigureWindowStyle", "docked")
%%
% Example Usage:
q0 = 1;  % Initial position
p0 = 0;  % Initial momentum (velocity)
dt = 0.01; % Time step
steps = 1000; % Number of integration steps
F_ext = @(q, p, t) 10;  % No external force

[q, p] = symplectic_spring_mass_damper(q0, p0, dt, steps, F_ext);


%%
function [q, p] = symplectic_spring_mass_damper(q0, p0, dt, steps, F_ext)
    % SYMPLECTIC_SPRING_MASS_DAMPER: 4th-order integration of damped mass-spring system
    % q0, p0: Initial conditions
    % dt: Time step
    % steps: Number of steps to integrate
    % F_ext: External force function handle (e.g., @(q,p,t) 0 for no force)

    % System parameters
    k = 1;  % Spring constant
    c = 1;  % Damping coefficient

    % 4th-order symplectic coefficients
    c_coeff = [1, 1] / (2 * (2 - 2^(1/3)));
    c_coeff = [c_coeff(1), (1 - 2^(1/3)) / (2 * (2 - 2^(1/3))), ...
               (1 - 2^(1/3)) / (2 * (2 - 2^(1/3))), c_coeff(2)];
    d_coeff = [1 / (2 - 2^(1/3)), -2^(1/3) / (2 - 2^(1/3)), ...
               1 / (2 - 2^(1/3)), 0];

    % Initialize
    q = zeros(1, steps);
    p = zeros(1, steps);
    q(1) = q0;
    p(1) = p0;

    % Time integration loop
    for step = 1:steps-1
        for i = 1:length(d_coeff)
            % Position update
            tempQ(step+1) = q(step) + c_coeff(i) * dt * p(step);
            
            % Momentum update (includes damping and external force)
            damping_force = -c * p(step);  
            spring_force = -k * q(step);
            external_force = F_ext(q, p, step * dt);
            tempP(step+1) = p(step) + d_coeff(i) * dt * (spring_force + damping_force + external_force);
        end
    end
end