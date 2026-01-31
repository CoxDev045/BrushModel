clear; close all; clc;

% --- Physical Parameters ---
k = 1e-4; c = 1e-7; m = 1e-10;
num_elements = 4;

% --- Time Stepping ---
omega_n = sqrt(k/m);
fs = 100 * (omega_n / (2*pi)); % 20x natural frequency for stability
dt = 1/fs;
time = 0:dt:0.5; % Simulation time
N = length(time);

% --- Newmark-beta coefficients (Average Acceleration) ---
gamma = 0.5; 
beta = 0.25;

% Precompute integration constants
a0 = 1/(beta*dt^2);
a1 = gamma/(beta*dt);
a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;
a5 = dt/2 * (gamma/beta - 2);

% Effective Stiffness (Since all elements are same, this is a scalar)
k_eff = k + a0*m + a1*c;

% --- Initial Conditions (Rows = States, Cols = Elements) ---
u = zeros(num_elements, N); % Displacement
v = zeros(num_elements, N); % Velocity
a = zeros(num_elements, N); % Acceleration

% Set initial displacement for each element (example)
u(:, 1) = [0.1, 0.2, 0.3, 0.4]; 

% Initial acceleration from F = ma -> a = (F - cv - ku)/m
% For simplicity, assume v0 = 0 and F0 = 0
a(:, 1) = (-c*v(:, 1) - k*u(:, 1)) / m;

% --- Forcing Functions ---
P = @(t) 0.0005 * sin(4 * t); % External drive
mu = @(vel) 0.5 * (tanh(vel-2)+1); 
 
% --- Integration Loop ---
for i = 1:N-1
    t_next = time(i+1);
    
    % 1. Calculate Forcing for each element
    % tau depends on current displacement (u)
    tau = mu(v(:, i)) * P(t_next); 
    u_max = tau / k;
    % 2. Effective Load (R_eff)
    % This is the core of Newmark-beta
    R_eff = tau + m*(a0*u(:, i) + a2*v(:, i) + a3*a(:, i)) ...
                + c*(a1*u(:, i) + a4*v(:, i) + a5*a(:, i));
    
    % 3. Solve for next displacement
    u_temp = R_eff / k_eff;
    % if u_temp >= u_max
    %     u_temp = u_max;
    % end
    u(:, i+1) = u_temp;   
    
    % 4. Update velocity and acceleration
    a_next = a0*(u(:, i+1) - u(:, i)) - a2*v(:, i) - a3*a(:, i);
    v_next = v(:, i) + (1-gamma)*dt*a(:, i) + gamma*dt*a_next;
    
    a(:, i+1) = a_next;
    v(:, i+1) = v_next;
end

% --- Visualization ---
T = tiledlayout(1,2);

nexttile
plot(time, u);
title('Displacement of 4 Independent Elements (Newmark-\beta)');
xlabel('Time (s)'); ylabel('Disp (\delta)');
legend('Elem 1','Elem 2','Elem 3','Elem 4');

nexttile
plot(u.', v.');
title('State Space of 4 Independent Elements (Newmark-\beta)');
xlabel('Disp (\delta)'); ylabel('Vel (\dot{\delta})')
legend('Elem 1','Elem 2','Elem 3','Elem 4');
%%
% --- Updated Parameters for 2D ---
k_vec = [1e-4; 1.2e-4]; % [kx; ky] - can be different
c_vec = [1e-7; 1.1e-7]; % [cx; cy]
m = 1e-10;
num_elements = 4;

% Newmark coefficients (remain the same as previous)
gamma = 0.5; 
beta = 0.25;

% Precompute integration constants
a0 = 1/(beta*dt^2);
a1 = gamma/(beta*dt);
a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;
a5 = dt/2 * (gamma/beta - 2);

% Effective Stiffness Matrix (2x2 for each element)
% Since X and Y are decoupled, this is a diagonal matrix
K_eff = diag([k_vec(1) + a0*m + a1*c_vec(1), ...
              k_vec(2) + a0*m + a1*c_vec(2)]);

% --- Initial Conditions (2 DOFs x num_elements x time) ---
u = zeros(2, num_elements, N); 
v = zeros(2, num_elements, N);
a = zeros(2, num_elements, N);

% Initial acceleration calculation (vectorized for 2 DOFs)
% a = (F - C*v - K*u) / m
for e = 1:num_elements
    u(:, e, 1) = [0.1; 0.1]; % Initial [dx; dy]
    % a = (Force - [cx 0; 0 cy]*v - [kx 0; 0 ky]*u) / m
    a(:, e, 1) = (-diag(c_vec)*v(:, e, 1) - diag(k_vec)*u(:, e, 1)) / m;
end

% --- Integration Loop ---
for i = 1:N-1
    % Forcing: tau_x and tau_y
    % You can define these as 2x1 vectors for each element
    P_t = [sin(100*time(i+1)); cos(80*time(i+1))];
    
    for e = 1:num_elements
        % mu depends on the specific displacement component
        mu_val = [mu(u(1, e, i)); mu(u(2, e, i))];
        tau = mu_val .* P_t; % Element-wise multiplication
        
        % Effective Load Vector (2x1)
        R_eff = tau + m*(a0*u(:, e, i) + a2*v(:, e, i) + a3*a(:, e, i)) ...
                    + diag(c_vec)*(a1*u(:, e, i) + a4*v(:, e, i) + a5*a(:, e, i));
        
        % Solve for next displacement (2x1)
        % Using \ (backslash) is efficient for this 2x2
        u(:, e, i+1) = K_eff \ R_eff;
        
        % Update v and a (2x1)
        a_next = a0*(u(:, e, i+1) - u(:, e, i)) - a2*v(:, e, i) - a3*a(:, e, i);
        v(:, e, i+1) = v(:, e, i) + (1-gamma)*dt*a(:, e, i) + gamma*dt*a_next;
        a(:, e, i+1) = a_next;
    end
end

% --- Visualization ---
T = tiledlayout(1,2);

nexttile
plot(time, squeeze(u(1, 1, :)));
hold on
plot(time, squeeze(u(1, 2, :)));
title('Displacement of 4 Independent Elements (Newmark-\beta)');
xlabel('Time (s)'); ylabel('Disp (\delta)');
% legend('Elem 1','Elem 2','Elem 3','Elem 4');
legend('X', 'Y')

nexttile
plot(squeeze(u(1, 1, :)), squeeze(u(2,1,:)));
title('State Space of 4 Independent Elements (Newmark-\beta)');
xlabel('Disp (\delta_x)'); ylabel('Disp (\delta_y)')
% legend('Elem 1','Elem 2','Elem 3','Elem 4');