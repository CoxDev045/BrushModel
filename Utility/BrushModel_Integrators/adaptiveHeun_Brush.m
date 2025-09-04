function [y_next, h, updated_obj] = adaptiveHeun_Brush(func, dt, t, X_vec, brush_obj)
    % ADAPTIVEHEUN_BRUSH performs one step of numerical integration.
    % An adaptive timestep with euler as first step and trapezoidal as
    % finer resolution step.
    %
    %   f: Function handle for the ODE: dy/dt = f(t, y)
    %   dt: Current time step
    %   t: Current time
    %   X_vec: Current solution vector
    %   args: Arguments for function func(t, x, args)
    %
    %   y_next: Solution at time t+h
    %   h: Time step at the end of the current step 
    y = X_vec;
    h = dt;
    hasToRunAgain = true;
    rtol = 1e-6;
    atol = 1e-6;
    h_min = 100 * eps;
    h_max = 1;
    max_iters = 5;

    % Calculate tolerance value
    tolerance = max(rtol * norm(y), atol);
    
    % Initialise updated_obj with brush_obj
    updated_obj = brush_obj;

    % Start counter
    iters = 0;
    while hasToRunAgain
        % --- Stage 1: RK1 Rule (Euler Step) ---
        % Solve for y_n+1 at t_n+1 = t + h
        [k1, euler_obj] = func(t, y, brush_obj);
        
        y_euler = X_vec + h * k1;
    
        % --- Stage 2: Trapezoidal rule (Heun's Method) ---
        % Solve for y_next at t_next = t + h
        [k2, trapz_obj] = func(t + h, y_euler, euler_obj);
        
        y_heun = y + 0.5 * (k1 + h * k2);
        % Estimate error betweem two schemes
        error = norm(y_heun - y_euler);
 
        % Calculate scaling factor
        S = 0.9 * (tolerance / error)^(1/2);
        
        if error <= tolerance
            % Accept current timestep
            y_next = y_heun;
            % Increase the step size
            h = min(h_max, S * h);
            hasToRunAgain = false;
            % Update brush_obj with new values
            updated_obj = trapz_obj;
        else
            % Reject the step
            h = max(h_min, S * h);
        end
        % Increment counter
        iters = iters + 1;
        % break out of loop to avoid infinite while loop
        % Can break out of loop if h == h_min because no refinement will
        % happen
        if iters >= max_iters || h == h_min
            hasToRunAgain = false;
            % Assume trapezoidal rule to be more accurate
            y_next = y_heun;
        end
    end

end