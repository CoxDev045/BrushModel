function [x_next, Fty_next, exitFlag, output] = solveTrustRegionDogLeg(func, t, x0, options)
    %SOLVETRUSTREGIONDOGLEG An implementation of the popular
    %trust-region-dogleg algorithm. This method makes use of a gradient
    %descent step and a newton step. It compares the two and chooses the
    %one that falls within a certain "region of trust", a circle with
    %radius rk that we know the solver did not overshoot the minimum.
    %Should the both solutions  fall outside the trust region, we solve the
    %quadratic equation to find the closest point on the trust-region
    %boundary and choose this as our next point. Function continues until
    %some exit condition is met.
    %
    % INPUTS
    % func: function handle of the function we want to minimise
    % t:    current value of t. Function is assumed to be a function of X
    %       and t
    % X0:   Initial guess at time t
    % options:  Struct containing the exit conditions of the solver:
    %           MaxIters:   Max iterations allowed
    %           JacTol:     Minimum value of the norm of the Jacobian
    %                       needed for us to know solution converged to a stationary
    %                       point
    %           FuncTol:    Minimum value of function evaluated at (t, X)
    %
    % OUTPUTS
    % X_next:   Solution returned by the solver. Point that minimises func
    % Fty_next: Function value at func(t, X)
    % exitFlag: Flag indicating exit condition
    % output:   Struct with additional parameters of interest:
    %           X_opt               Value that minimises the supplied function
    %           ExitFlag            Flag indication exit condition
    %           Jacobian            Jacobian at optimised point
    %           Hessian             Hessian at optimised point
    %           FunVal              Function value at optimised point
    %           TrustRegionRadius   Final radius of trust-region
    %           ImprovementRatio    Final improvement ratio (Actual improvement / Quadratic model of function)
    %           Iterations          Iterations took to optimise
    %           message             Exit message
    % 
    %----------------------------------------------------------------------
    %                       Function Starts here
    %----------------------------------------------------------------------
    exitFlag = 0;
    % Initialise point xk with x0
    xk = x0;
    % Choose initial radius rk
    rk = 1;
    rk_max = 2;
    rk_min = eps;
    % Choose threshold value η ∈ [0, 0.25)
    eta = 0.0625;
    
    % Extract options
    thresh_scal = options.ThreshScal;
    % Evaluate function at inital value
    Fty = func(t, xk);
    
    % Calculate Jacobian at point xk
    [J,fac] = numjac(func,t,xk,Fty,thresh_scal);
    % Calculate Hessian at point xk
    Hk = J.' * J;
    % Calculate norm of Jac
    normJac = norm(J);
    iters = 0;
    while exitFlag == 0
        % % % Compute steepest descent step
        % % p_sd = - J.' * Fty;
        % % 
        % % if rcond(Hk) >= 1e-3
        % %     % Compute Gauss-Newton step
        % %     p_gn = decomposition(Hk) \ p_sd;
        % % else
        % %     p_gn = pinv(Hk) * p_sd;
        % % end

        % Calculate gradient
        g = J.' * Fty;

        % Check conditioning of Hessian
        Hk_cond = rcond(Hk);
        
        if Hk_cond >= 1e-3
            % Compute unconstrained minimiser step
            p_gn = -1 * decomposition(Hk) \ g;
        
        elseif Hk_cond < 1e-6 % Apply a pre-conditioner
            % Get diagonal. This will serve as the pre-conditioner
            D = diag(Hk);
            D = D + 1e-3;
            
            % Create pre-conditioned Hessian
            Hk_prec = Hk ./ D;

            % Solve pre-conditioned system
            g_prec = g ./ D;
            p_gn_prec = -decomposition(Hk_prec) \ g_prec;
            p_gn = p_gn_prec;
        else
            p_gn = -1 * pinv(Hk) * g;
        end
        % Compute steepest descent step
        p_sd = -1 * (g.'* g) / (g.' * Hk * g) * g; 

        % Compute norm of steps
        norm_p_sd = norm(p_sd);
        norm_p_gn = norm(p_gn);

        if norm_p_gn <= rk 
            % Gauss-Newton step is within or on the trust region
            pk = p_gn;
        elseif norm_p_sd >= rk
            % Steepest descent an Gauss-Newton is outside trust-region
            % Scale steepest descent point to be on trust region circle
            pk = (rk / norm_p_sd) * p_sd;
        else
            % The step is on the dogleg path
            tau = solveQuadraticProblem(p_sd, p_gn, rk); % Find tau such that the combined step is on the boundary
            pk = p_sd + (tau - 1) * (p_gn - p_sd);
        end
        % Calculate new point
        x_next = xk + pk;
        % Evaluate function at new point
        Fty_next =func(t, x_next);
         % Calculate actual reduction
        ActualRed = 0.5 * norm(Fty)^2 - 0.5 * norm(Fty_next)^2;
        % Calculate predicted reduction
        PredRed =  - pk.' * g - 0.5 * pk.' * Hk * pk;
        if PredRed <=0
            % Avoid division by zero or negative reduction
            rho_k = 0;
        else
            % Compute reduction ratio
            rho_k = ActualRed / PredRed;
        end
        % Check if step is smaller than some threshold eta
        if rho_k > eta
            % If larger, accept step
            xk = x_next;
            Fty = Fty_next;
            % Calculate Jacobian at new point
            [J_next,fac_next] = numjac(func,t,x_next,Fty_next,thresh_scal, fac);   
            % Update Jacobian
            J = J_next;
            fac = fac_next;
            % Calculate exit condition
            normJac = norm(J);
            % Update Hessian approximation
            Hk = J.' * J;
        end

        if rho_k < 0.25
            % Model is bad, shrink radius
            rk = max(rk_min, 0.25 * rk);
        elseif rho_k > 0.75 && norm(pk) == rk
            % Model is good, grow radius
            rk = min(2 * rk, rk_max);
        end
        
        iters = iters +1;

        % Check exit conditions
        [exitFlag, msg] = checkExit(normJac, rho_k, iters, Fty, options);
        
    end
    output = struct('X_opt', x_next,...
                    'ExitFlag', exitFlag,...
                    'Jacobian', J,...
                    'Hessian', Hk,...
                    'FunVal', Fty,...
                    'TrustRegionRadius', rk,...
                    'ImprovementRatio', rho_k,...
                    'Iterations', iters,...
                    'message', msg);

    function [exitFlag, msg] = checkExit(normJac, rho_k, iters, FunVal, options)
        exitFlag = 0;
        msg = '';
        if rho_k <= 0
            exitFlag = -2;
            msg = sprintf('Solver could not reduce the function any further! Local minimum possible.\n');
        elseif iters >= options.MaxIters
            exitFlag = -1;
            msg = sprintf('Failed to converge within max iterations.\n');
        elseif normJac <= options.JacTol
            exitFlag = 1; 
            msg = sprintf('Successfully reduced Jacobian to within tolerance of %g\n', options.JacTol);
        elseif FunVal <= options.FunTol
            exitFlag = 2;
            msg = sprintf('Successfully reduced function value to within tolerance of %g\n', options.FunTol);
        end
        
    end

    function tau = solveQuadraticProblem(p_sd, p_gn, rk)

        a = (p_gn.' * p_gn) - 2 * (p_gn.' * p_sd) + (p_sd.' * p_sd);
        b = 2 * (p_gn.' * p_sd) - 2 * (p_sd.' * p_sd);
        c = (p_sd.' * p_sd) - rk^2;
        
        top_pos = -b + sqrt( (b)^2 - 4 * (a) * (c) );
        bot = 2 * a;
    
        tau = 1 + top_pos / bot;
        
        % Handle potential floating-point errors
        if imag(tau) ~= 0 || tau < 1 || tau > 2
            % No valid solution found, should not happen with correct inputs
            error('No valid value of tau could be found!');
        end
    end
end

