function X_next = evaluateYoshida(func, dt, t, X_vec)
    % 4th Order Yoshida Integrator for systems with explicit acceleration function.
    % This function performs one integration step from t to t+dt.
    %
    % Inputs:
    %   func:       Function handle to compute acceleration.
    %               Expected signature: a = func(t, x, v, args_for_func{:})
    %               where x is position and v is velocity (potentially half-step velocity).
    %   dt:         Time step.
    %   t:          Current time.
    %   X_vec:      Current state vector [x_current; v_current; a_current].
    %               x_current: Position(s) at time t.
    %               v_current: Velocity(s) at time t.
    %               a_current: Acceleration(s) at time t (calculated from previous step).
    %   varargin:   Optional additional arguments to pass directly to func.
    %
    % Output:
    %   X_next:     Next state vector [x_new; v_new; a_new].
    %               x_new: Position(s) at time t+dt.
    %               v_new: Velocity(s) at time t+dt.
    %               a_new: Acceleration(s) at time t+dt.

    % --- 1. Define Correct Yoshida Coefficients ---
    % These are for the common A-B-A composition (position-velocity-position)
    % which requires 4 position updates and 3 velocity updates.
    p = 2^(1/3);
    w0 = -p / (2 - p); % approx -0.7657
    w1 = 1 / (2 - p);  % approx 1.3512

    c1 = w1 / 2; c4 = c1;
    c2 = (w0 + w1) / 2; c3 = c2;
    
    d1 = w1; d3 = d1;
    d2 = w0;

    % --- 2. Robust State Variable Extraction ---
    % Assuming X_vec = [x_components; v_components; a_old_components]
    N = numel(X_vec); % Use numel for total elements
    num_dims = N / 3; % Number of dimensions (e.g., 1 for 1D, 2 for 2D)

    if mod(N, 3) ~= 0 || N < 3
        error('Input state vector X_vec must have a total number of elements divisible by 3 (position, velocity, and acceleration components for each dimension). For 1D, N=3; for 2D, N=6, etc.');
    end
    
    x = X_vec(1:num_dims);            % Position(s) at time t
    v = X_vec(num_dims+1 : 2*num_dims); % Velocity(s) at time t
    a_old = X_vec(2*num_dims+1 : 3*num_dims); % Acceleration(s) at time t

    % --- 3. Integration Steps (Composition) ---
    % Here, 'A' step is x update, 'B' step is v update, 'func' provides 'a' (acceleration)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c1 * dt)
    x = x + c1 .* v .* dt;
    
    % B(d1 * dt) - calculate 'a' at new position (x) and current velocity (v)
    % The 'v' here (v_current) is often used for damping in the 'a' calculation
    a_temp1 = get_acceleration(func,t, [x; v]); % Calculate a(t+c1*dt) approximately
    v = v + d1 .* a_temp1 .* dt; % Update velocity

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c2 * dt)
    x = x + c2 .* v .* dt;

    % B(d2 * dt)
    a_temp2 = get_acceleration(func,t, [x; v]); % Calculate a(t+(c1+c2)*dt) approximately
    v = v + d2 .* a_temp2 .* dt; % Update velocity

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c3 * dt)
    x = x + c3 .* v .* dt;

    % B(d3 * dt)
    a_temp3 = get_acceleration(func,t, [x; v]); % Calculate a(t+(c1+c2+c3)*dt) approximately
    v = v + d3 .* a_temp3 .* dt; % Update velocity

    %%%%%%%%%%%%%%%%%%%%%%%%% Step 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A(c4 * dt) - Final position update
    x = x + c4 .* v .* dt;

    % --- 4. Recalculate Final Acceleration for Output Consistency ---
    % After all position and velocity updates, calculate the acceleration
    % corresponding to the final (t+dt) state.
    a_final = get_acceleration(func,t+dt, [x; v]); % a(t+dt)

    % --- 5. Prepare Output State Vector ---
    % X_next contains [new_position; new_velocity; new_acceleration]
    X_next = [x; v; a_final];
end

function acceleration = get_acceleration(func, t_in, X_vec)
    % This function adapts the state-space dynamics to return only acceleration.
    % It's designed to be passed as 'func' to integrate_VelocityVerlet_general.
    % Call the original state-space dynamics function
    dX = func(t_in, X_vec);
    
    [N, ~] = size(dX);

    if mod(N, 2) == 0
        if N > 2
            % Extract the acceleration component (dvx is the second element of dX)
            acceleration = dX(3:end); 
            acceleration = acceleration(:);
        else
            acceleration = dX(2);
        end
    else
        acceleration = [];
        error('Output contains wrong number of variables! Please ensure that output contains at least one velocity and one acceleration component!');
    end
end