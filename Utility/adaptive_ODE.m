function [y_next, h_next] = adaptive_ODE(func, dt, t_current, t_target, X_vec, args)
    % adaptive_ODE performs one step of numerical integration.
    % An adaptive timestep with RK4 as first step and RK5 as
    % finer resolution step.
    %
    %   f: Function handle for the ODE: dy/dt = f(t, y)
    %   dt: Current time step
    %   t_current: Current time
    %   t_target: Final time step to solve for
    %   X_vec: Current solution vector
    %   args: Arguments for function func(t, x, args)
    %
    %   y_next: Solution at time t+h
    %   h_next: Time step at the end of the current step 

    % --- RKF45 Coefficients ---
    % Coefficients for the k_i stages
    A = [0; 1/2; 1/2; 1; 2/3; 2/10];
    C = [0     , 0       , 0      , 0     , 0;
         1/2   , 0       , 0      , 0     , 0;
         1/4   , 1/4     , 0      , 0     , 0;
         0     , -1      , 2      , 0     , 0;
         7/27  , 10/27   , 0      , 1/27  , 0;
         28/625, -125/625, 546/625, 54/625, -378/625];
    
    % Coefficients for the 4th and 5th order solutions
    B_4 = [1/6, 0, 4/6, 3, 0, 0];
    B_5 = [14/336, 0, 0, 35/336, 162/336, 125/336];

    % Initialise h and y
    y = X_vec;
    h_current = dt;
    
    % Adaptive time step parameters
    rtol = 1e-6;
    atol = 1e-6;
    h_min = 100 * eps;
    h_max = 0.1 * abs(t_target);
    max_retries = 5;
    tolerance = max(rtol * norm(y), atol);

    retries = 0;
    while retries < max_retries
        % Calculate all k_i values
        k1 = h_current * func(t_current + A(1)*h_current, y, args{:});
        k2 = h_current * func(t_current + A(2)*h_current, y + C(2,1)*k1, args{:});
        k3 = h_current * func(t_current + A(3)*h_current, y + C(3,1)*k1 + C(3,2)*k2, args{:});
        k4 = h_current * func(t_current + A(4)*h_current, y + C(4,1)*k1 + C(4,2)*k2 + C(4,3)*k3, args{:});
        k5 = h_current * func(t_current + A(5)*h_current, y + C(5,1)*k1 + C(5,2)*k2 + C(5,3)*k3 + C(5,4)*k4, args{:});
        k6 = h_current * func(t_current + A(6)*h_current, y + C(6,1)*k1 + C(6,2)*k2 + C(6,3)*k3 + C(6,4)*k4 + C(6,5)*k5, args{:});
        
        % Calculate the two solutions (4th and 5th order)
        y_RK4 = y + B_4(1)*k1 + B_4(2)*k2 + B_4(3)*k3 + B_4(4)*k4 + B_4(5)*k5 + B_4(6)*k6;
        y_RK5 = y + B_5(1)*k1 + B_5(2)*k2 + B_5(3)*k3 + B_5(4)*k4 + B_5(5)*k5 + B_5(6)*k6;
        
        % Estimate the local error
        error = norm(y_RK5 - y_RK4);
        
        % Calculate tolerance and scaling factor
        S = 0.9 * (tolerance / error)^(1/5);
        
        if error <= tolerance
            % Step Accepted
            y_next = y_RK5; % Use the higher-order solution
            h_next = min(h_max, h_current * S);
            return; % Exit the function after a successful step.
        else
            % Step Rejected
            h_current = max(h_min, h_current * S);
            retries = retries + 1;
        end
    end
    
    % If max_retries is reached, return last solution with a warning
    warning('Adaptive step failed to converge within max retries.');
    y_next = y; % Return the original point
    h_next = h_current;
end
