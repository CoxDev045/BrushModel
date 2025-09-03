function X_next = evaluateImplicitEuler(func, dt, t, X_vec)
    % evaluateImplicitEuler_Newton: Implements the Implicit Euler method using MATLAB's fsolve
    % (which uses a trust-region or quasi-Newton algorithm).
    % This approach is more robust for stiff ODEs compared to fixed-point iteration.
    %
    % Inputs:
    %   func:       Function handle for the ODE. Signature: dXdt = func(time, state_vector, args{:})
    %   dt:         Time step size.
    %   t:          Current time (t_n).
    %   X_vec:      Current state vector (X_n).
    %   varargin:   Optional additional arguments to pass to func.
    %
    % Output:
    %   X_next:     Approximation of the state vector at time t + dt (X_{n+1}).

    % Ensure 'optimoptions' and 'fsolve' are available in the MATLAB environment.
    % This code requires the Optimization Toolbox.

    t_next = t + dt;

    % --- Define the function to find the root of ---
    % We want to solve for X_next such that:
    % X_next - X_vec - dt * func(t_next, X_next, varargin{:}) = 0
    % Define this as an anonymous function for fsolve.
    
    equation_to_solve = @(t, X_guess) norm( X_guess - X_vec - dt * func(t_next, X_guess) );

    % --- Initial Guess for X_next ---
    % A good initial guess is crucial for fsolve convergence.
    % The Explicit Euler step is often used as a starting point.
    initial_guess = X_vec + dt * func(t, X_vec);

    % --- Set fsolve options (optional, but recommended for control) ---
    % You can adjust these based on your specific problem's needs.
    % 'Display': 'off', 'final', 'iter'
    % 'FunctionTolerance': Tolerance on the function value
    % 'StepTolerance': Tolerance on the change in X_next
    % options = optimoptions('fsolve', ...
    %                        'Display', 'off', ... % Suppress verbose output from fsolve
    %                        'FunctionTolerance', 1e-8, ... % Tolerance for F(X_next) close to zero
    %                        'StepTolerance', 1e-8);     % Tolerance for change in X_next
    % options = optimset('Display', 'off', ... % Suppress verbose output from fsolve
    %                    'TolFun', 1e-8,...
    %                    'TolX', 1e-8); % Tolerance for change in X_next

    options.JacTol = 1e-8;
    options.ThreshScal = sqrt(eps) * ones(size(X_vec));
    options.MaxIters = 1000;
    options.FunTol = 1e-8;

    % --- Call fsolve to solve for X_next ---
    % fsolve returns X_next, and potentially fval (value of the function at solution),
    % exitflag, output structure, and Jacobian. We only need X_next for this function.
    % [X_next, ~, exitflag, output] = fminsearch(equation_to_solve, initial_guess, options);
    [X_next, ~, exitflag, output] = solveTrustRegionDogLeg(equation_to_solve, t, initial_guess, options);

    % % --- Check fsolve exit flag (optional, but good practice) ---
    % if exitflag <= 0 % exitflag < 1 typically means no convergence
    %     warning('evaluateImplicitEuler_Newton:fsolveNoConvergence', ...
    %             'fsolve did not converge successfully for t=%f, dt=%f. Exit flag: %d. Message: %s', ...
    %             t, dt, exitflag, output.message);
    % end
end