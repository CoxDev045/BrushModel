function [X_next, updated_obj] = integrateDynamics(func, dt, t, X_vec, method_name, obj, args)
%INTEGRATEDYNAMICS is a simple wrapper function that takes in the user's
%pre-defined dynamics as a function handle and integrates it
%forward in time using the method specified in the method_name input
    % INPUTS
    % func:         Function handle of user's dynamics. Assumed to be
    %               only a function of (t, X).
    % dt:           Time step used. Will return next value at t + dt
    % t:            Current time value
    % X_vec:        Current state vector
    % method_name:  String containing the method of choice
    %
    % OUTPUTS
    % X_next:       The integrated states at time (t + dt)
    %

    arguments
        func            
        dt              (1,1) single
        t               (1,1) single
        X_vec           (:,1) single
        method_name     
        obj
        args
    end

    if isa(func, 'function_handle')
        switch lower(method_name)
            case 'euler'
                [X_next, updated_obj] = evaluateEuler_Brush(func, dt, t, X_vec, obj, args);
            case 'implicit_euler'
                X_next = evaluateImplicitEuler(func, dt, t, X_vec);
            case 'adaptive_heun'
                t_current = t;
                t_target = t + dt;
                h_current = dt;
                while t_current < t_target
                    % Calculate required step
                    h_current = min(h_current,  t_target - t_current);
            
                    % Call the adaptive step function
                    [X_next, h_next] = adaptiveHeun(func, h_current, t_current, X_vec);
            
                    % Update current time based on the step taken
                    t_current = t_current + h_current;
                    % Update time step
                    h_current = h_next;
                    % Update solution
                    X_vec = X_next;
                end
            case 'verlet'
                X_next = evaluateVerlet(func, dt, t, X_vec);
            case 'velocity_verlet'
                X_next = evaluateVelocityVerlet(func, dt, t, X_vec);
            case 'tr_bdf2'
                X_next = evaluateTRBDF2(func, dt, t, X_vec);
            case 'ode23s'
                X_sol = ode23s(func, dt, t, X_vec);
                X_next = X_sol(end, :);
            case 'ode23tb'
                X_sol = ode23tb(func, dt, t, X_vec);
                X_next = X_sol(end, :);
            case 'ode23t'
                X_sol = ode23t(func, dt, t, X_vec);
                X_next = X_sol(end, :);
            case 'rk4'
                [X_next, updated_obj] = evaluateRK4_Brush(func, dt, t, X_vec, obj, args);
            case 'rkf5'
                X_next = evaluateRKF5(func, dt, t, X_vec);
            case 'adaptive_rk45'
                t_current = t;
                t_target = t + dt;
                h_current = dt;
                while t_current < t_target
                    % Calculate required step
                    h_current = min(h_current,  t_target - t_current);
            
                    % Call the adaptive step function
                    [X_next, h_next] = adaptiveRK45(func, h_current, t_current, X_vec);
            
                    % Update current time based on the step taken
                    t_current = t_current + h_current;
                    % Update time step
                    h_current = h_next;
                    % Update solution
                    X_vec = X_next;
                end
            case 'ode45'
                X_sol = ode45(func, dt, t, X_vec);
                X_next = X_sol(end, :);
            otherwise
                error('Unrecognised integration method: %s', method_name);
        end
    else
        error('Unrecognised function:  %s. Please provide a valid function handle!', func2str(func));
    end
end
