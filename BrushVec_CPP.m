classdef BrushVec_CPP %#codegen -args
    % BrushVec - A vectorized brush dynamics and friction model
    % Optimised to work with codegeneration to CPP using MATLAB
    % 
    % This class implements a physics-based model for brush dynamics with friction,
    % supporting multiple brush elements simultaneously through vectorization.
    %
    % The model handles both adhesion (sticking) and sliding states, with appropriate
    % physics for each regime. Various numerical integration methods are available.
    %
    % Example:
    %   brush = BrushVec(xPositions, yPositions, pressure, slidingVelocity);
    %   brush = brush.update_brush(omega, omega_z, re, v0, alpha, dt);
    
    properties (SetAccess = public)
        % Brush Model properties (Static)
        % % phi         (1,1) double = 1;%0.32;      % Anisotropy coefficient
        kx          (1,1) single = 1e-5;%8.85431878608975;            % Base x-stiffness
        ky          (1,1) single = 1e-4;%0.405027674510023;%0.37;      % Base y-stiffness
        cx          (1,1) single = 1e-5;%6.89967576251185;%1.78e-7;   % x-damping coefficient
        cy          (1,1) single = 1e-4;%0.00704256221258847;%1.40e-4;   % y-damping coefficient
        m_inv       (1,1) single = 936.605399758573;%1.3089005e+09;  % Inverse of Mass property
        m           (1,1) single = 0.00106768549514851;%7.64e-10;  % Mass

        % Friction Model Properties (Static)
        mu_0        (1,1) single = 0.00;                             % Static friction coefficient
        mu_m        (1,1) single = 1.20;      % Maximum friction coefficient
        h           (1,1) single = 1.5;%0.4;      % Friction model parameter
        p_0         (1,1) single = 0.02;      % Minimum pressure threshold
        p_ref       (1,1) single = 0.247645523118594;%0.39;      % Reference pressure
        q           (1,1) single = 0.390845735345209;%0.28;      % Pressure exponent
        v_m         (1,1) single = 5;%30.206780784527050;%23.47;     % Reference velocity

        % Dynamic properties (Changing throughout simulation)
        numBrushes  (1, 1) uint16          % Number of brushes in model
        x           (400,1) single         % x-coordinate of the brush
        minX        (1, 1)  single         % Minimum x-value that x can be in when brush is in contact patch
        maxX        (1, 1)  single         % Maximum x_value that x can be when brush is in contact patch
        y           (400,1) single         % y-coordinate of the brush
        delta_x     (400,1) single         % x-displacement of the brush
        delta_y     (400,1) single         % y-displacement of the brush
        tauX        (400,1) single         % x shear force
        tauY        (400,1) single         % y shear force
        mu          (400,1) single      % Friction coefficient
        press       (400,1) single         % Pressure
        vrx         (400,1) single         % Relative x velocity
        vry         (400,1) single         % Relative y velocity
        vs          (400,1) single         % Sliding velocity magnitude
        vs_x        (400,1) single         % X Sliding Velocity
        vs_y        (400,1) single         % Y Sliding Velocity
        prev3_vrx   (400,1) single         % Relative x velocity, 3 steps back
        prev3_vry   (400,1) single         % Relative y velocity, 3 steps back
        prev2_vrx   (400,1) single         % Relative x velocity, 2 steps back
        prev2_vry   (400,1) single         % Relative y velocity, 2 steps back
        prev1_vrx   (400,1) single         % Relative x velocity, 1 step back
        prev1_vry   (400,1) single         % Relative y velocity, 1 step back
        dvrx        (400,1) single         % Approximated Gradient of relative velocity in x direction
        dvry        (400,1) single         % Approximated Gradient of relative velocity in y direction
        theta_2     (400,1) single         % Angle between sliding velocity and horizontal
        theta_1     (400,1) single         % Angle between shear stress and horizontal
        slide       (400,1) logical        % Sliding state (true/false)
        passed      (400,1) logical        % To check whether the brush has moved past contact length

        % Integration method state variables
        ax          (400,1) single         % x-acceleration
        ay          (400,1) single         % y-acceleration
        ax_new      (400,1) single         % New x-acceleration for Verlet
        ay_new      (400,1) single         % New y-acceleration for Verlet
        prev1_vx    (400,1) single         % X Deformation Velocity 1 step back
        prev1_vy    (400,1) single         % Y Deformation Velocity 1 step back
        vx          (400,1) single         % x-velocity
        vy          (400,1) single         % y-velocity

        % --- RKF4(5) - III Method ---
        % Coefficients for the k_i stages
        A = [0; 1/4; 3/8; 12/13; 1; 1/2];
        C = [0        , 0         , 0         , 0        , 0;
             1/4      , 0         , 0         , 0        , 0;
             3/32     , 9/32      , 0         , 0        , 0;
             1932/2197, -7200/2197, 7296/2197 , 0        , 0;
             439/216  , -8        , 3680/513  , -845/4104, 0;
             -8/27    , 2         , -3544/2565, 1859/4104, -11/40];
            
        % Coefficients for the 4th and 5th order solutions
        B_4 = [25/216, 0, 1408/2565 , 2197/4104  , -1/5, 0];
        B_5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    end
    
    methods
        % Constructor
        function [obj] = BrushVec_CPP(xVal, yVal, press, numBrushes)
            % Constructor for BrushVec class
            %
            % Parameters:
            %   xVal  - x-coordinates of brush elements
            %   yVal  - y-coordinates of brush elements
            %   press - Vertical pressure at each element
            %   vs    - Sliding velocity at each element
            arguments (Input)
                xVal        (:,1)  single
                yVal        (:,1)  single
                press       (:,1)  single
                numBrushes  (1, 1) uint16  % **Must be a constant** 
            end
            
            % Grid size
            obj.numBrushes = numBrushes;

            % Validate inputs
            assert(size(xVal, 1) <= numBrushes);
            assert(size(yVal, 1) <= numBrushes);
            assert(size(press, 1) <= numBrushes);

            % Find these parameters seperately from data
            obj.h = 2.0;     % Stribeck Exponent
            obj.q = 0.390845735345209;     % Pressure Saturation Exponent
            obj.v_m = 30.206780784527050;   % Stribeck Velocity
            
            % Initialize properties
            obj.x = xVal;
            obj.minX = min(xVal);
            obj.maxX = max(xVal);
            obj.y = yVal;
            obj.press = press;                       % Vertical Pressure

        end


        function [obj] = solve_brush_slide_ode(obj, omega, omega_z, re, v0, alpha, dt)
            % Computes brush dynamics in sliding region
            %
            % Parameters:
            %   omega  - Angular velocity
            %   omega_z - Angular velocity around z-axis
            %   re     - Effective radius
            %   v0     - Initial velocity
            %   alpha  - Angle
            %   dt     - Time step
            %   integrationMethod - Method for numerical integration (optional)
            
            arguments
                obj
                omega   (:,1) single
                omega_z (1,1) single
                re      (1,1) single
                v0      (:,1) single
                alpha   (1,1) single
                dt      (1,1) single
            end
            
            assert(size(omega, 1) <= 400);
            assert(size(v0, 1) <= 400);

            % Update velocity history
            obj.prev3_vrx  = obj.prev2_vrx ;
            obj.prev3_vry  = obj.prev2_vry ;
            obj.prev2_vrx  = obj.prev1_vrx ;
            obj.prev2_vry  = obj.prev1_vry ;
            obj.prev1_vrx  = obj.vrx ;
            obj.prev1_vry  = obj.vry ;

            % Get Sliding index
            slideInd = obj.slide;
            % Construct current solution vectors for use with integration
            X_vec = [obj.delta_x(slideInd); obj.vx(slideInd)];
            Y_vec = [obj.delta_y(slideInd); obj.vy(slideInd)];
            
            % Apply selected integration method
            stillRunning = true;
            while stillRunning
                % Integrate in X-direction until target time is reached
                if tx_current < t_target
                    [X_next, hx_next] = integrateODE(slidingDynamicsX, X_vec, h_current, t_current);
                    tx_current = tx_current + hx_next;
                end
                % Integrate in Y-direction until target time is reached
                if ty_current < t_target 
                    [Y_next, hy_next] = integrateODE(slidingDynamicsY, Y_vec, h_current, t_current);
                    ty_current = ty_current + hy_next;
                end
                % Check if both directions have reached the target time
                if tx_current >= t_target && ty_current >= t_target
                    stillRunning = false;
                end
            end
            obj.delta_x(slideInd) = X_next(1:obj.numBrushes, 1);
            obj.vx(slideInd)      = X_next(obj.numBrushes+1:end, 1);
            obj.delta_y(slideInd) = Y_next(obj.numBrushes, 1);
            obj.vy(slideInd)      = Y_next(obj.numBrushes+1:end, 1); 
           
        end

        function [obj] = calculateStressesAdaptive()
            obj.vrx(slideInd) = omega .* re + omega_z .* (obj.y(slideInd) + obj.delta_y(slideInd)) - v0 .* cos(alpha);
            obj.vry(slideInd) = -omega_z .* (obj.x(slideInd) + obj.delta_x(slideInd)) - v0 .* sin(alpha);
            
            % Calculate sliding velocity
            obj.vs_y(slideInd) = real( obj.vy(slideInd) - obj.vry(slideInd) );
            obj.vs_x(slideInd) = real( obj.vx(slideInd) - obj.vrx(slideInd) );
            obj.vs(slideInd) = real( hypot(obj.vs_x(slideInd), obj.vs_y(slideInd)) ) .* sign(obj.vs_x(slideInd));

            obj.theta_1(slideInd) = real( atan2(obj.vs_y(slideInd), obj.vs_x(slideInd)) ); 
            obj.theta_2(slideInd) = obj.theta_1(slideInd) - pi;

            % % cos_theta_1 = cos(obj.theta_1(slideInd));
            % % mask = abs(cos_theta_1) > 1e-8;
            % % obj.vs(slideInd) = mask .* (obj.vx(slideInd) - obj.vrx(slideInd)) ./ cos_theta_1;

            % Calculate friction coefficient
            obj.mu = obj.compute_friction_coefficient(obj.press, obj.p_0, obj.p_ref, obj.q, obj.mu_0, obj.mu_m, obj.h, obj.vs, obj.v_m);

            % Calculate stresses
            obj.tauX(slideInd) = obj.mu(slideInd) .* obj.press(slideInd) .* cos(obj.theta_2(slideInd));
            obj.tauY(slideInd) = obj.mu(slideInd) .* obj.press(slideInd) .* sin(obj.theta_2(slideInd));

        end

        function dX = slidingDynamicsX(t,X, SlideInd)
            obj.dDeltax(SlideInd) = obj.vx(SlideInd);
            obj.dVx(SlideInd) = (obj.tauX(SlideInd) - obj.kx(SlideInd) .* obj.delta_x(SlideInd) - obj.cx(SlideInd) .* obj.vx(SlideInd)) ./ obj.m;
            dX = [obj.dDeltax(SlideInd);
                  obj.dVx(SlideInd)];            
        end

        function dY = slidingDynamicsY(t, Y)

            obj.dDeltay(SlideInd) = obj.vy;
            obj.dVy(SlideInd) = (obj.tauY(SlideInd) - obj.ky(SlideInd) .* obj.delta_y(SlideInd) - obj.cy(SlideInd) .* obj.vy(SlideInd)) ./ obj.m;
            dY = [obj.dDeltay(SlideInd);
                  obj.dVy(SlideInd)];
        end

        function [y_next, h_next] = integrateODE(func, X0, dt, t_current)
            % adaptive_ODE performs one step of numerical integration.
            % An adaptive timestep with RK4 as first step and RK5 as
            % finer resolution step.
            %
            %   func: Function handle for the ODE: dy/dt = f(t, y)
            %   X_vec: Current solution vector
            %   dt: Current time step
            %   t_current: Current time
            %
            %   y_next: Solution at time t+h
            %   h_next: Time step at the end of the current step 

            % Initialise h and X0
            X_vec = X0;
            h_current = dt;
            
            % Adaptive time step parameters
            rtol = 1e-14;
            atol = 1e-13;
            h_min = 100 * eps;
            h_max = 1e-3;
            max_retries = 20;
            tolerance = max(rtol * norm(X_vec), atol);
        
            retries = 0;
            while retries < max_retries
                % Calculate all k_i values
                k1 = h_current * func(t_current + obj.A(1)*h_current, X_vec);
                k2 = h_current * func(t_current + obj.A(2)*h_current, X_vec + obj.C(2,1)*k1);
                k3 = h_current * func(t_current + obj.A(3)*h_current, X_vec + obj.C(3,1)*k1 + obj.C(3,2)*k2);
                k4 = h_current * func(t_current + obj.A(4)*h_current, X_vec + obj.C(4,1)*k1 + obj.C(4,2)*k2 + obj.C(4,3)*k3);
                k5 = h_current * func(t_current + obj.A(5)*h_current, X_vec + obj.C(5,1)*k1 + obj.C(5,2)*k2 + obj.C(5,3)*k3 + obj.C(5,4)*k4);
                k6 = h_current * func(t_current + obj.A(6)*h_current, X_vec + obj.C(6,1)*k1 + obj.C(6,2)*k2 + obj.C(6,3)*k3 + obj.C(6,4)*k4 + obj.C(6,5)*k5);
                % Calculate the two solutions (4th and 5th order)
                y_RK4 = X_vec + obj.B_4(1)*k1 + obj.B_4(2)*k2 + obj.B_4(3)*k3 + obj.B_4(4)*k4 + obj.B_4(5)*k5 + obj.B_4(6)*k6;
                y_RK5 = X_vec + obj.B_5(1)*k1 + obj.B_5(2)*k2 + obj.B_5(3)*k3 + obj.B_5(4)*k4 + obj.B_5(5)*k5 + obj.B_5(6)*k6;
      
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
            y_next = X_vec; 
            h_next = h_current;

        end

        function [obj] = integrate_verlet(obj, dt)
            % Corrected Velocity Verlet integration method
            %
            % This function uses an object-oriented approach where state is
            % stored as properties of the object 'obj'.
        
            % Get Sliding index
            slideInd = obj.slide;
            
            % Pre-calculate common expressions
            half_dt = 0.5 .* dt;
            
            % --- 1. Store old acceleration ---
            % obj.ax will be overwritten, so we need to save its value for the
            % second half-step velocity update.
            ax_old = obj.ax(slideInd);
            ay_old = obj.ay(slideInd);
            
            % --- 2. First half-step velocity update ---
            % v(t + dt/2) = v(t) + a(t) * (dt/2)
            % Note: Your code used obj.prev1_vx here. I've used obj.vx as it
            % is the current velocity v(t) being updated.
            v_half_vx = obj.vx(slideInd) + ax_old .* half_dt;
            v_half_vy = obj.vy(slideInd) + ay_old .* half_dt;
        
            % --- 3. Full position update ---
            % x(t + dt) = x(t) + v(t + dt/2) * dt
            % This is the key change: REMOVE the a*dt^2/2 term
            obj.delta_x(slideInd) = obj.delta_x(slideInd) + v_half_vx .* dt;
            obj.delta_y(slideInd) = obj.delta_y(slideInd) + v_half_vy .* dt;
            
            % --- 4. Calculate new accelerations (at time t+dt) ---
            % a(t + dt) = Force(x(t+dt), v(t+dt/2)) / m
            % We use the NEW position and the HALF-STEP velocity for the force calculation.
            ax_verlet = (obj.tauX(slideInd) - obj.kx .* obj.delta_x(slideInd) - obj.cx .* v_half_vx) * obj.m_inv;
            ay_verlet = (obj.tauY(slideInd) - obj.ky .* obj.delta_y(slideInd) - obj.cy .* v_half_vy) * obj.m_inv;
            
            % --- 5. Second half-step velocity update with new accelerations ---
            % v(t + dt) = v(t + dt/2) + a(t + dt) * (dt/2)
            % This is the standard, cleaner Velocity Verlet formula.
            obj.vx(slideInd) = v_half_vx + ax_verlet .* half_dt;
            obj.vy(slideInd) = v_half_vy + ay_verlet .* half_dt;
            
            % --- 6. Update acceleration state for the next step ---
            % Store the new acceleration so it can be used as the 'old' acceleration in the next step.
            obj.ax(slideInd) = ax_verlet;
            obj.ay(slideInd) = ay_verlet;
        end
        
        function [obj] = solve_stick_ode(obj, omega, omega_z, re, v0, alpha, dt)
            % Computes brush dynamics in adhesion/sticking region
            %
            % Parameters:
            %   omega  - Angular velocity
            %   omega_z - Angular velocity around z-axis
            %   re     - Effective radius
            %   v0     - Initial velocity
            %   alpha  - Angle
            %   dt     - Time step
            
            arguments
                obj
                omega   (:,1) single
                omega_z (1,1) single
                re      (1,1) single
                v0      (:,1) single
                alpha   (1,1) single
                dt      (1,1) single
            end
            
            % Update velocity history
            obj.prev3_vrx   = obj.prev2_vrx  ;
            obj.prev3_vry   = obj.prev2_vry  ;
            obj.prev2_vrx   = obj.prev1_vrx  ;
            obj.prev2_vry   = obj.prev1_vry  ;
            obj.prev1_vrx   = obj.vrx  ;
            obj.prev1_vry   = obj.vry  ;

            % Calculate New Velocities
            obj.vrx(~obj.slide) = omega .* re + omega_z .* (obj.y(~obj.slide) + obj.delta_y(~obj.slide)) - v0 .* cos(alpha);
            obj.vry(~obj.slide) = -omega_z .* (obj.x(~obj.slide) + obj.delta_x(~obj.slide)) - v0 .* sin(alpha);

            % Update displacements for sticking
            obj.delta_x(~obj.slide) = obj.delta_x(~obj.slide) + obj.vrx(~obj.slide) .* dt;
            obj.delta_y(~obj.slide) = obj.delta_y(~obj.slide) + obj.vry(~obj.slide) .* dt;

            % 2nd order backward difference for acceleration estimation
            obj.dvrx(~obj.slide) = (3 .* obj.vrx(~obj.slide) - 4 .* obj.prev1_vrx(~obj.slide) + obj.prev2_vrx(~obj.slide)) ./ (2 * dt);
            obj.dvry(~obj.slide) = (3 .* obj.vry(~obj.slide) - 4 .* obj.prev1_vry(~obj.slide) + obj.prev2_vry(~obj.slide)) ./ (2 * dt);

            % Calculate shear stress
            obj.tauX(~obj.slide) = obj.kx .* obj.delta_x(~obj.slide) + obj.cx .* obj.vrx(~obj.slide) + obj.m .* obj.dvrx(~obj.slide);
            obj.tauY(~obj.slide) = obj.ky .* obj.delta_y(~obj.slide) + obj.cy .* obj.vry(~obj.slide) + obj.m .* obj.dvry(~obj.slide);
        end

        function [obj] = update_brush(obj, pressVal, omega, omega_z, re, v0, alpha, dt)
            % Updates brush state for one time step
            %
            % Parameters:
            %   omega  - Angular velocity
            %   omega_z - Angular velocity around z-axis
            %   re     - Effective radius
            %   v0     - Initial velocity
            %   alpha  - Angle
            %   dt     - Time step
            %   integrationMethod - Numerical integration method (optional)
            arguments
                obj
                pressVal    (400,1) single
                omega       (1,1)   single
                omega_z     (1,1)   single
                re          (1,1)   single
                v0          (1,1)   single
                alpha       (1,1)   single
                dt          (1,1)   single
            end
            
            % Validate inputs
            validateattributes(omega, {'single'}, {'scalar'});
            validateattributes(omega_z, {'single'}, {'scalar'});
            validateattributes(re, {'single'}, {'scalar'});
            validateattributes(v0, {'single'}, {'scalar'});
            validateattributes(alpha, {'single'}, {'scalar'});
            validateattributes(dt, {'single'}, {'scalar', 'positive'});
            
            % Assume uniform angular velocity and translational velocity of
            % road element is constant across contact patch
            omega_vec = repmat(omega, size(pressVal));
            v0_vec = repmat(v0, size(pressVal));

            %%%%%%%%%%%%%% Update Pressure Dependent Properties %%%%%%%%%%%%%%%%
            obj.press = pressVal;
            %%%%%%%%%%%%% Remove Pressure Dependent stiffness for now %%%%%
            % % obj.ky = obj.ky_0 + obj.ky_l .* obj.press;
            % % obj.kx = obj.phi .* obj.ky;
            
            % Determine if elements are sliding or sticking
            obj.slide = obj.is_sliding(obj.tauX, obj.tauY, obj.mu_0, obj.press);
            
            % Update dynamics based on state
            if any(obj.slide)
                % Handle partial sliding
                if ~all(obj.slide)
                    % Split object to handle slide and stick separately
                    
                    if ~isempty(omega_vec(obj.slide))
                        slide_obj = obj;
                        % Process sliding elements
                        slide_obj = slide_obj.solve_brush_slide_ode(omega_vec(obj.slide), omega_z, re, v0_vec(obj.slide), alpha, dt);
                        % Combine results back
                        obj.delta_x(obj.slide) = slide_obj.delta_x(obj.slide);
                        obj.delta_y(obj.slide) = slide_obj.delta_y(obj.slide);
                        obj.vx(obj.slide) = slide_obj.vx(obj.slide);
                        obj.vy(obj.slide) = slide_obj.vy(obj.slide);

                        obj.vrx(obj.slide) = slide_obj.vrx(obj.slide);
                        obj.vry(obj.slide) = slide_obj.vry(obj.slide);

                        obj.tauX(obj.slide) = slide_obj.tauX(obj.slide);
                        obj.tauY(obj.slide) = slide_obj.tauY(obj.slide);
                    end

                    if ~isempty(omega_vec(~obj.slide))
                        stick_obj = obj;
                        % Process sticking elements
                        stick_obj = stick_obj.solve_stick_ode(omega_vec(~obj.slide), omega_z, re, v0_vec(~obj.slide), alpha, dt);
                        % Combine results back
                        obj.delta_x(~obj.slide) = stick_obj.delta_x(~obj.slide);
                        obj.delta_y(~obj.slide) = stick_obj.delta_y(~obj.slide);

                        obj.vrx(~obj.slide) = stick_obj.vrx(~obj.slide);
                        obj.vry(~obj.slide) = stick_obj.vry(~obj.slide);

                        obj.tauX(~obj.slide) = stick_obj.tauX(~obj.slide);
                        obj.tauY(~obj.slide) = stick_obj.tauY(~obj.slide);
                    end
                else
                    % All elements are sliding
                    obj = obj.solve_brush_slide_ode(omega_vec, omega_z, re, v0_vec, alpha, dt);
                end
            else
                % All elements are sticking
                obj = obj.solve_stick_ode(omega_vec, omega_z, re, v0_vec, alpha, dt);
            end
            
            obj = obj.updateXvalues(omega, re, dt);
            
        end
    
        function [obj] = updateXvalues(obj, omega, re, dt)
            arguments
                obj
                omega   (1,1) single
                re      (1,1) single
                dt      (1,1) single
            end
            % Update x position
            obj.x = obj.x - omega .* re .* dt;
    
            % check if brush has passed through contact patch
            hasPassed = obj.passedThrough(obj.x, obj.minX);
            obj.passed = hasPassed;
            if any(hasPassed)
                % fprintf('Brush has passed when omega = %.6f \n', omega);
                % Set the properties of the brushes that has passed
                % back to zero
                %%% Bring brush back to front edge %%%
                obj.x(hasPassed) = obj.maxX;
                %%% Reset Displacements %%%
                obj.delta_x(hasPassed) = 0;
                obj.delta_y(hasPassed) = 0;
                %%% Reset Stresses %%%%
                obj.tauX(hasPassed) = obj.mu(hasPassed) .* obj.press(hasPassed) .* cos(-pi);
                obj.tauY(hasPassed) = 0;
                %%% Reset friction coefficient
                obj.mu(hasPassed) = obj.mu(hasPassed);
                %%% Reset sliding velocity
                obj.vs(hasPassed) = 0;
                %%% Reset sliding condition %%%
                obj.slide(hasPassed) = true;
                %%% Reset deformation velocities %%%
                obj.vx(hasPassed) = 0;
                obj.vy(hasPassed) = 0;
                %%% Reset relative velocities %%%
                obj.vrx(hasPassed) = 0;
                obj.vry(hasPassed) = 0;

                obj.prev3_vrx(hasPassed) = 0;          % Relative x velocity, 3 steps back
                obj.prev3_vry(hasPassed) = 0;         % Relative y velocity, 3 steps back
                obj.prev2_vrx(hasPassed) = 0;         % Relative x velocity, 2 steps back
                obj.prev2_vry(hasPassed) = 0;         % Relative y velocity, 2 steps back
                obj.prev1_vrx(hasPassed) = 0;         % Relative x velocity, 1 step back
                obj.prev1_vry(hasPassed) = 0;         % Relative y velocity, 1 step back
                obj.theta_2(hasPassed) = 0;           % Angle between sliding velocity and horizontal
                obj.theta_1(hasPassed) = 0;           % Angle between shear stress and horizontal

                % Integration method state variables
                obj.ax(hasPassed) = 0;         % x-acceleration
                obj.ay(hasPassed) = 0;         % y-acceleration
                obj.ax_new(hasPassed) = 0;     % New x-acceleration for Verlet
                obj.ay_new(hasPassed) = 0;     % New y-acceleration for Verlet
            end
        end
    end


    methods (Static)
        function [slide] = is_sliding(tauX, tauY, mu_0, press)
            % Determines if elements are in sliding state
            %
            % Parameters:
            %   tauX  - X-component of shear stress
            %   tauY  - Y-component of shear stress
            %   mu_0  - Static friction coefficient
            %   press - Pressure
            %
            % Returns:
            %   slide - Logical array indicating sliding elements
            
            arguments
                tauX   (400,1)  single
                tauY   (400,1)  single
                mu_0   (1,1)    single
                press  (400,1)  single
            end
            
            tau = hypot(tauX, tauY);
            slide = (tau >= mu_0 .* press);
        end
        
        function [mu] = compute_friction_coefficient(press, p_0, p_ref, q, mu_0, mu_m, h, vs, v_m)
            % Computes the dynamic friction coefficient
            %
            % Parameters:
            %   press - Pressure
            %   p_0   - Minimum pressure threshold
            %   p_ref - Reference pressure
            %   q     - Pressure exponent
            %   mu_0  - Static friction coefficient
            %   mu_m  - Maximum friction coefficient
            %   h     - Friction model parameter
            %   vs    - Sliding velocity
            %   v_m   - Reference velocity
            %
            % Returns:
            %   mu    - Dynamic friction coefficient
            
            arguments
                press   (400,1) single
                p_0     (1,1)   single
                p_ref   (1,1)   single
                q       (1,1)   single
                mu_0    (1,1)   single
                mu_m    (1,1)   single
                h       (1,1)   single
                vs      (400,1) single
                v_m     (1,1)   single
            end
            
            % Calculate pressure-dependent term
            masked_p = max(press, p_0);
            sat_p = real( (masked_p ./ p_ref) .^(-q) );

            % Ensure velocity is valid for logarithm
            vm_safe = max(abs(v_m), eps) .* sign(v_m);
            vs_safe = max(abs(vs), eps) .* sign(vs);
            
            % Calculate velocity-dependent term
            logterm = log10(abs(vs_safe ./ vm_safe));
            scaleTerm = (mu_m - mu_0);
            expTerm = exp(-h.^2 .* logterm.^2);
            mu = sat_p .* (mu_0 +  scaleTerm .* expTerm) .* sign(vs_safe);

            % Check for NaNs in pressure-related calculations
            if any(isnan(sat_p))
                warning('NaN detected in sat_p! Check press, p_ref, p_0, and q.');
            end
            
            % Check for NaNs in velocity-related calculations
            if any(isnan(logterm))
                warning('NaN detected in logterm! Check vs and v_m.');
            end
            
            % Check final friction coefficient calculation
            if any(isnan(mu))
                warning('NaN detected in mu! Check sat_p, scaleTerm, and expTerm.');
            end

        end

        function [passed] = passedThrough(currentX, minX)
            % Determines if elements have passed through the contact patch
            %
            % Parameters:
            %   currentX  - current X value
            %   minX  - minimum value X can be when brush is in contact patch
            %
            % Returns:
            %   passed - Logical array indicating passed elements
            arguments
                currentX    (400,1) single
                minX        (400,1) single
            end
            passed = currentX <= minX;
        end

    end
end

            