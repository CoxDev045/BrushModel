function [y_next] = trbdf2_step(func, dt, t_val, X_vec, args)
    % TRBDF2_STEP performs one step of the TR-BDF2 integration scheme.
    %
    %   f: Function handle for the ODE: dy/dt = f(t, y)
    %   t: Current time
    %   y: Current solution vector
    %   h: Timestep size
    %   tol: Tolerance for the Newton's method
    %   max_iter: Maximum iterations for Newton's method
    %
    %   y_next: Solution at time t+h
    %   t_next: Time at the end of the step (t+h)
    
    % The optimal parameter gamma for L-stability
    gamma = 2 - sqrt(2);
    
    % --- Stage 1: Trapezoidal Rule (Implicit Step) ---
    % Solve for y_intermediate at t_gamma = t + gamma*h
    % The equation to solve is: y_intermediate - y - (gamma*h/2)*(f(t,y) + f(t_gamma, y_intermediate)) = 0
    y = X_vec;
    y_intermediate = X_vec; % Initial guess for Newton's method
    h = dt;
    % Calculate intermediate timestep
    t_intermediate = t_val + gamma * h;
    
    F1 = @(t, y_i_gamma) y_i_gamma - ( y + (gamma * h / 2) * ( func(t_val, y, args{:}) + func(t_intermediate, y_i_gamma, args{:}) ) );

    % % % --- Set fsolve options (optional, but recommended for control) ---
    % % % You can adjust these based on your specific problem's needs.
    % % % 'Display': 'off', 'final', 'iter'
    % % % 'FunctionTolerance': Tolerance on the function value
    % % % 'StepTolerance': Tolerance on the change in X_next
    % % % options = optimoptions('lsqnonlin', ...
    % % %                        'Display', 'off', ... % Suppress verbose output from fsolve
    % % %                        'FunctionTolerance', 1e-8, ... % Tolerance for F(X_next) close to zero
    % % %                        'StepTolerance', 1e-8);     % Tolerance for change in X_next
    options.JacTol = 1e-3;
    options.ThreshScal = sqrt(eps) * ones(size(y));
    options.MaxIters = 500;
    options.FunTol = 1e-8;

    % % --- Call fsolve to solve for X_next ---
    % % fsolve returns X_next, and potentially fval (value of the function at solution),
    % % exitflag, output structure, and Jacobian. We only need X_next for this function.
    % [y_intermediate, ~, exitflag, output] = fsolve(F1, y_intermediate, options);
    [y_intermediate, ~, exitflag, output] = solveTrustRegionDogLeg(F1, t_val, y_intermediate, options);
    % [y_intermediate, ~, ~, exitflag, output] = lsqnonlin(F1, y_intermediate, [], [], options);

    % --- Check fsolve exit flag (optional, but good practice) ---
    if exitflag <= 0 % exitflag < 1 typically means no convergence
        warning('trbdf2_step:solverNoConvergence', ...
                'Solver did not converge successfully for t=%f, dt=%f. Exit flag: %d. Message: %s', ...
                t_val, dt, exitflag, output.message);
    end
    
    % --- Stage 2: BDF2 (Implicit Step) ---
    % Solve for y_next at t_next = t + h
    % The equation to solve is: y_next - (1/(2-gamma))*( ((1-gamma)^2/gamma)*y + (1/gamma)*y_intermediate) - (1-gamma)/(2-gamma)*h*f(t_next, y_next) = 0
    
    y_next = y_intermediate; % Initial guess for Newton's method
    F2 = @(t, y_ip1) y_ip1 - (1/(2-gamma)) * ( (1/gamma)*y_intermediate - ((1-gamma)^2/gamma)*y + (1-gamma)*h*func(t_val + h, y_ip1, args{:}) );
   
    % %  % --- Call fsolve to solve for X_next ---
    % % % fsolve returns X_next, and potentially fval (value of the function at solution),
    % % % exitflag, output structure, and Jacobian. We only need X_next for this function.
    % % [y_next, ~, exitflag, output] = fsolve(F2, y_next, options);
    [y_next, ~, exitflag, output] = solveTrustRegionDogLeg(F2, t_val, y_next, options);
    % [y_next, ~, ~, exitflag, output] = lsqnonlin(F2, y_next, [], [], options);

    % --- Check fsolve exit flag (optional, but good practice) ---
    if exitflag <= 0 % exitflag < 1 typically means no convergence
        warning('trbdf2_step:solverNoConvergence', ...
                'Solver did not converge successfully for t=%f, dt=%f. Exit flag: %d. Message: %s', ...
                t_val, dt, exitflag, output.message);
    end

end