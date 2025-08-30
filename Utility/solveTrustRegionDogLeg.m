function [x_next, Fty_next, exitFlag, output] = solveTrustRegionDogLeg(func, t, x0, options)
    %SOLVETRUSTREGIONDOGLEG Summary of this function goes here
    exitFlag = 0;
    % Initialise point xk with x0
    xk = x0;
    % Choose initial radius rk
    rk = 1;
    rk_max = 10;
    rk_min = 1e-8;
    % Choose threshold value η ∈ [0, 0.25)
    eta = 0.125;
    
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
        % Compute steepest descent step
        p_sd = - J.' * Fty;
        % Compute Gauss-Newton step
        p_gn = decomposition(Hk) \ p_sd;

        % Compute norm of steps
        norm_p_sd = norm(p_sd);
        norm_p_gn = norm(p_gn);

        if norm_p_gn <= rk 
            % Gauss-Newton step is within or on the trust region
            pk = p_gn;
        elseif norm_p_sd >= rk
            % Steepest descent is outside trust-region
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
        ActualRed = norm(Fty)^2 - norm(Fty_next)^2;
        % Calculate predicted reduction
        PredRed = -(J.' * Fty ).' * pk - 0.5 * pk.' * Hk * pk;
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
        [exitFlag, msg] = checkExit(normJac, iters, Fty, options);
        
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

    % fprintf(msg);

end

function tau = solveQuadraticProblem(p_sd, p_gn, rk)
    a_vec = p_gn - p_sd;
    b_vec = p_sd;

    a = a_vec.' * a_vec;
    b = 2 * (a_vec.' * b_vec);
    c = (b_vec.' * b_vec - rk^2);
    
    top_pos = -b + sqrt( (b)^2 - 4 * (a) * (c) );
    bot = 2 * a;

    tau = 1 + top_pos / bot;
    
    % Handle potential floating-point errors
    if imag(tau) ~= 0 || tau < 1 || tau > 2
        % No valid solution found, should not happen with correct inputs
        error('No valid value of tau could be found!');
    end
    
end

function [exitFlag, msg] = checkExit(normJac, iters, FunVal, options)
    exitFlag = 0;
    msg = '';
    if iters >= options.MaxIters
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