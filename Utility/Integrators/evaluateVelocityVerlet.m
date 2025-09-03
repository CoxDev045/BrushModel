function X_next = evaluateVelocityVerlet(func, dt, t, X_vec)
%EVALUATEVELOCITYVERLET
% Generalized Velocity Verlet integration method.
% This function computes one step of the integration.
%
    % Inputs:
    %   accel_func: Function handle to compute acceleration.
    %               Expected signature: a = accel_func(t, x, v, args_for_accel_func{:})
    %               where x is position and v is velocity (potentially half-step velocity).
    %   dt:         Time step.
    %   t:          Current time.
    %   X_vec:      Current state vector [x_current; v_current; a_old].
    %               x_current is position at time t.
    %               v_current is velocity at time t.
    %               a_old is acceleration at time t (from previous step's calculation).
    %   varargin:   Optional additional arguments to pass directly to accel_func.
    %
    % Output:
    %   X_next:     Next state vector [x_new; v_new; a_new].
    %               x_new is position at time t+dt.
    %               v_new is velocity at time t+dt.
    %               a_new is acceleration at time t+dt.

    % Extract current state variables
    [N, ~] = size(X_vec);
    if mod(N, 3) == 0
        if N == 6
            % Extract the acceleration component (dvx is the second element of dX)
            x_current = X_vec(1:2); 
            x_current = x_current(:);
            v_current = X_vec(3:4); 
            v_current = v_current(:);
            a_old = X_vec(5:6); 
            a_old = a_old(:);
            
        else
            x_current = X_vec(1); % x(t)
            v_current = X_vec(2); % v(t)
            a_old     = X_vec(3); % a(t) - acceleration from the previous step

        end
    else
        x_current = []; % x(t)
        v_current = []; % v(t)
        a_old     = []; % a(t) - acceleration from the previous step

        error('Output contains wrong number of variables! Please ensure that output contains at least one velocity and one acceleration component!');
    end
    % 1. Calculate velocity at half-step
    v_half_step = v_current + a_old * (dt / 2);

    % 2. Update position
    x_new = x_current + v_half_step * dt;

    % 3. Calculate new acceleration using the provided accel_func
    % Pass the new position (x_new) and the half-step velocity (v_half_step)
    % along with any additional parameters to the acceleration function.
    X_vec_for_func = [x_new; v_half_step];
    a_new = get_acceleration(func, t, X_vec_for_func);% a(t + dt

    % 4. Calculate full-step velocity
    v_new = v_half_step + a_new * (dt / 2); % v(t + dt)
    
    % X_next contains [new_position; new_velocity; new_acceleration]
    X_next = [x_new; v_new; a_new]; 
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
