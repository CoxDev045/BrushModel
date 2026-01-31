clear;close all;clc
k = 1e-7;
m = 1e-10;

omega = sqrt(k/m);
zeta = 0;

A = [0,1;
    -omega^2,-zeta * omega];

dt = 0.001;
M_FE = eye(size(A)) + dt * A;
M_BE = inv(eye(size(A)) - dt * A);
M_CN = (eye(size(A)) - dt/2 * A) \ (eye(size(A)) + dt/2 * A);
M_RK4 = eye(size(A)) + dt * A + ((dt * A)^2)/2 + ((dt * A)^3)/6 + ((dt * A)^4)/24;
M_VV = [(2-dt^2*omega^2)/2, (2*dt-dt^2*zeta*omega)/2;
        -dt*omega^2*(4-dt^2*omega^2)/(4*(1+dt*zeta*omega/2)), (2-dt*zeta*omega-dt^2*omega^2-dt^3*omega^3*zeta/2)/(2*(1+dt*zeta*omega/2))];

gamma = 2 - sqrt(2);
K1 = (eye(size(A)) - gamma * dt * A / 2);
K2 = (eye(size(A)) + gamma * dt * A / 2);
M1 = K1 \ K2;
c1 = 1/(gamma * (2-gamma));
c2 = (1-gamma)^2 / (gamma * (2-gamma));
M_TRBDF2 = K1 \ (c1 * M1 - c2 * eye(size(A)));

num_Ints = 6;

fprintf('Eigenvalues of update matrix: FE\n')
beta = eig(M_FE);
disp(beta)
fprintf('Magnitude of eigenvalues: %.16f\n', abs(beta(1)))

% fprintf('Eigenvalues of update matrix: BE\n')
% beta = eig(M_BE);
% disp(beta)
% 
% fprintf('Magnitude of eigenvalues: %.6g\n', abs(beta(1)))
% 
% fprintf('Eigenvalues of update matrix: CN\n')
% beta = eig(M_CN);
% disp(beta)
% 
% fprintf('Magnitude of eigenvalues: %.6g\n', abs(beta(1)))

t_final = 10;
t = linspace(0, t_final, t_final/dt);


X = zeros(2 * num_Ints, length(t));
vel_ind = linspace(1,num_Ints,num_Ints) * 2;
disp_ind = linspace(1,num_Ints,num_Ints) * 2 -1;

X(disp_ind, 1) = 2e0;
fk = 1e-3;
M_BE_inv = (eye(size(A)) - dt * A);
for i = 1:length(t)-1
    % Calculate forcing vector at current state
    force_vec = [0; -fk * sign(X(2, i))];
    X(1:2, i+1) = M_FE * X(1:2, i) + dt * force_vec;
    
    force_vec = [0; -fk * sign(X(4, i))];
    X(3:4, i+1) = M_BE_inv \ X(3:4, i) + dt * force_vec;
    
    force_vec = [0; -fk * sign(X(6, i))];
    X(5:6, i+1) = M_CN * X(5:6, i) + dt * force_vec;
    
    force_vec = [0; -fk * sign(X(8, i))];
    X(7:8, i+1) = M_RK4 * X(7:8, i) + dt * force_vec;
    
    force_vec = [0; -fk * sign(X(10, i))];
    X(9:10, i+1) = M_VV * X(9:10, i) + dt * force_vec;
    
    force_vec = [0; -fk * sign(X(12, i))];
    X(11:12, i+1) = M_TRBDF2 * X(11:12, i) + dt * force_vec;
end
F_ext = @forcing_func;
initial_state = [X(1, 1); X(2, 1)];
args = {omega, zeta, F_ext};

my_dynamics = @(t, X) springMassDamperDynamics(t, X, args{:});
options = odeset("RelTol",1e-6, "AbsTol",1e-6,"Stats","on"); 
tic;
[t_sol, X_sol] = ode23s(my_dynamics,t , initial_state, options);
fprintf('Elapsed time: %gs \n', toc)

%%
KE = 0.5 * m * X(vel_ind, :).^2;         % Kinetic Energy
PE = 0.5 * k * X(disp_ind, :).^2;    % Potential Energy (elastic)
TotalMechanicalEnergy = KE + PE;      % Total Mechanical Energy

T = tiledlayout(2,2);
T.Padding = "compact";
T.TileSpacing = "tight";

nexttile
hold on
for i = 1:num_Ints
    plot(t, X(disp_ind(i), :))
end
grid on
xlabel('Time [s]')
ylabel('Displacement [m]')
legend('FE', 'BE', 'CN', 'RK4', 'VV', 'TRBF2')

nexttile
hold on
for i = 1:num_Ints
    plot(t, X(vel_ind(i), :))
end
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('FE', 'BE', 'CN', 'RK4', 'VV', 'TRBDF2')

nexttile
plot(t, TotalMechanicalEnergy)
grid on
xlabel('Time [s]')
ylabel('Energy [J]')
title('Total Mechanical Energy')
legend('FE', 'BE', 'CN', 'RK4', 'VV', 'TRBDF2')

nexttile
plot(KE.', PE.')
grid on
xlabel('Time [s]')
ylabel('Energy [J]')
title('Total Mechanical Energy')
legend('FE', 'BE', 'CN', 'RK4', 'VV', 'TRBDF2')


%%
function dX = springMassDamperDynamics(t, X, omega, zeta, Forcing)

    dx = X(2, :);
    dvx = (Forcing(t, X) - omega^2 .* X(1, :) - zeta * omega .* X(2, :));
    dX = [dx; dvx];
end

function F = forcing_func(t, X)
    % This is your state-dependent forcing function.
    % It calculates the forces based on time and the current system state.
    v = X(2, :);

    % Example: Forcing that depends on velocity
    F_magnitude = 10; % Example magnitude
    F = F_magnitude * sign(v); % Example: a force that opposes velocity
end