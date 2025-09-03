function X_next = evaluateVerlet(func, dt, t, X_vec)
    % Position Verlet integration method for dynamics.
    % X_vec should be [current_position; previous_position]
   
    % Extract current state variables
    [N, ~] = size(X_vec);
    if mod(N, 3) == 0
        if N == 6
            % Extract the acceleration component (dvx is the second element of dX)
            x_current = X_vec(1:2); 
            x_current = x_current(:);
            v_estimated = X_vec(3:4); 
            v_estimated = v_estimated(:);
            x_previous = X_vec(5:6); 
            x_previous = x_previous(:);
            
        else
            x_current   = X_vec(1); % x(t)
            v_estimated = X_vec(2); % v(t)
            x_previous  = X_vec(3); % x(t-dt) - displacement from the previous step

        end
    else
        error('Output contains wrong number of variables! Please ensure that output contains at least one velocity and one acceleration component!');
    end

    % 1. Calculate acceleration at current position
    X_vec_for_func = [x_current;v_estimated];
    a_current = get_acceleration(func, t, X_vec_for_func); % a(t)

    % 2. Update position
    x_new = 2 * x_current - x_previous + a_current * dt.^2; % x(t + dt)

    % 3. Calculate velocity for output (not used for propagation in core Verlet)
    % A better estimate for v(t+dt) for analysis is often (x(t+dt) - x(t-dt)) / (2*dt)
    v_new_estimated = (x_new - x_previous) / (2 * dt); % v(t + dt)

    % X_next contains [new_position;estimated_new_velocity; current_position (for next step's "previous"); ]
    X_next = [x_new; v_new_estimated; x_current]; 
end