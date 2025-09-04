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
        kx          (1,1) single = 8.85431878608975;            % Base x-stiffness
        ky          (1,1) single = 0.405027674510023;%0.37;      % Base y-stiffness
        cx          (1,1) single = 6.89967576251185;%1.78e-7;   % x-damping coefficient
        cy          (1,1) single = 0.00704256221258847;%1.40e-4;   % y-damping coefficient
        m_inv       (1,1) single = 1;%936.605399758573;%1.3089005e+09;  % Inverse of Mass property
        m           (1,1) single = 1;%0.00106768549514851;%7.64e-10;  % Mass

        % Friction Model Properties (Static)
        mu_0        (1,1) single = 0.01;                        % Static friction coefficient
        mu_m        (1,1) single = 1.20;                        % Maximum friction coefficient
        h           (1,1) single = 0.75;%1.5;%0.4;                    % Friction model parameter
        p_0         (1,1) single = 0.02;                        % Minimum pressure threshold
        p_ref       (1,1) single = 0.247645523118594;%0.39;     % Reference pressure
        q           (1,1) single = 0.390845735345209;%0.28;     % Pressure exponent
        v_m         (1,1) single = 2.5;%5;%30.206780784527050;%23.47;% Reference velocity

        % Dynamic properties (Changing throughout simulation)
        numBrushes      (1, 1) uint16   = 400;       % Number of brushes in model
        x               (400,1) single  = 0;       % x-coordinate of the brush
        minX            (1, 1)  single  = -1;       % Minimum x-value that x can be in when brush is in contact patch
        maxX            (1, 1)  single  = 1;       % Maximum x_value that x can be when brush is in contact patch
        y               (400,1) single  = 0;       % y-coordinate of the brush
        delta_x         (400,1) single  = 0;       % x-displacement of the brush
        delta_y         (400,1) single  = 0;       % y-displacement of the brush
        tauX            (400,1) single  = 0;       % x shear force
        tauY            (400,1) single  = 0;       % y shear force
        mu              (400,1) single  = 0.01;       % Friction coefficient
        press           (400,1) single  = 0;       % Pressure
        vrx             (400,1) single  = 0;       % Relative x velocity
        vry             (400,1) single  = 0;       % Relative y velocity
        vs              (400,1) single  = 0;       % Sliding velocity magnitude
        vs_x            (400,1) single  = 0;       % X Sliding Velocity
        vs_y            (400,1) single  = 0;       % Y Sliding Velocity
        prev1_delta_x   (400,1) single  = 0;       % Relative x displacement, 1 step back
        prev1_delta_y   (400,1) single  = 0;       % Relative y displacement, 1 step back
        prev3_vrx       (400,1) single  = 0;       % Relative x velocity, 3 steps back
        prev3_vry       (400,1) single  = 0;       % Relative y velocity, 3 steps back
        prev2_vrx       (400,1) single  = 0;       % Relative x velocity, 2 steps back
        prev2_vry       (400,1) single  = 0;       % Relative y velocity, 2 steps back
        prev1_vrx       (400,1) single  = 0;       % Relative x velocity, 1 step back
        prev1_vry       (400,1) single  = 0;       % Relative y velocity, 1 step back
        dvrx            (400,1) single  = 0;       % Approximated Gradient of relative velocity in x direction
        dvry            (400,1) single  = 0;       % Approximated Gradient of relative velocity in y direction
        theta_2         (400,1) single  = 0;       % Angle between sliding velocity and horizontal
        theta_1         (400,1) single  = 0;       % Angle between shear stress and horizontal
        slide           (400,1) logical = false;       % Sliding state (true/false)
        passed          (400,1) logical = false;       % To check whether the brush has moved past contact length

        % Integration method state variables
        ax          (400,1) single = 0;        % x-acceleration
        ay          (400,1) single = 0;       % y-acceleration
        ax_new      (400,1) single = 0;        % New x-acceleration for Verlet
        ay_new      (400,1) single = 0;        % New y-acceleration for Verlet
        prev1_vx    (400,1) single = 0;        % X Deformation Velocity 1 step back
        prev1_vy    (400,1) single = 0;        % Y Deformation Velocity 1 step back
        vx          (400,1) single = 0;        % x-velocity
        vy          (400,1) single = 0;        % y-velocity

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
            obj.h = 0.75;%2.0;     % Stribeck Exponent
            obj.q = 0.390845735345209;     % Pressure Saturation Exponent
            obj.v_m = 2.5;%30.206780784527050;   % Stribeck Velocity

            
            % Initialize properties from inputs
            obj.x = xVal;
            obj.minX = min(xVal);
            obj.maxX = max(xVal);
            obj.y = yVal;
            obj.press = press;         % Vertical Pressure
            % Initialise other properties that are not zero
            obj.mu = obj.mu_0;
            

        end


        function [obj] = solve_brush_slide_ode(obj, omega, omega_z, re, v0, alpha, dt, t)
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
                t       (1,1) single
            end
            
            assert(size(omega, 1) <= 400);
            assert(size(v0, 1) <= 400);

            % Get Sliding index
            slideInd = obj.slide;
            % Construct current solution vectors for use with integration
            X_vec = [obj.delta_x(slideInd);
                     obj.delta_y(slideInd);
                     obj.vx(slideInd);
                     obj.vy(slideInd)];
            forcing_args = struct('omega', omega,... 
                                  're', re,... 
                                  'omega_z', omega_z,... 
                                  'v0', v0,... 
                                  'alpha', alpha);

            % brush_dynamics = @ BrushVec_CPP.slidingDynamics(t,X, obj, forcing_args);
            
            [X_next, slide_obj] = integrateDynamics(@slidingDynamics, dt, t, X_vec, 'euler', obj, forcing_args);
            
            %
            num_masses = length(X_next) / 4;
            % Extract into object
            obj.delta_x(slideInd) = X_next(1:num_masses, 1);
            obj.vx(slideInd)      = X_next(num_masses+1:2 * num_masses, 1);
            obj.delta_y(slideInd) = X_next(2 * num_masses+1: 3 * num_masses, 1);
            obj.vy(slideInd)      = X_next(3 * num_masses+1:4 * num_masses, 1); 

            % Load the other states from slide_obj to current obj
            obj.mu(slideInd) = slide_obj.mu(slideInd);
            obj.vs(slideInd) = slide_obj.vs(slideInd);
            obj.vs_x(slideInd) = slide_obj.vx(slideInd);
            obj.vs_y(slideInd) = slide_obj.vy(slideInd);
           
        end

        function [obj] = solve_stick_ode(obj, omega, omega_z, re, v0, alpha, dt)
            %SOLVE_STICK_ODE is a method that steps the simulation forward
            %by one time step. This method is only called when there are
            %brushes that are not sliding. This means there is only
            %relative acceleration between road element and brush due to
            %the relative velocity. The stresses are calculated from the
            %new displacements, new velocities and new accelerations.
            % (d^δ/dt^2 = d^2v_r/dt
            % => 𝜏 = k * δ + c * v_r + m * dv_r/dt)
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

            forcing_args = struct( 'omega', omega,...  
                                   'omega_z', omega_z,...
                                   're', re, ...
                                   'v0', v0, ...
                                   'alpha', alpha, ...
                                   'dt', dt);

            % brush_dynamics = @(t, X) BrushVec_CPP.stickingDynamics(t, X, obj, forcing_args);

            % Get Sliding index
            slideInd = ~obj.slide;
            % Construct current solution vectors for use with integration
            X_vec = [obj.delta_x(slideInd);
                     obj.delta_y(slideInd);
                     obj.vrx(slideInd);
                     obj.vry(slideInd)];

            X_next = integrateDynamics(@stickingDynamics, dt, single(0), X_vec, 'euler', obj, forcing_args);

            num_masses = length(X_next) / 4;
            % Extract into object
            obj.delta_x(~obj.slide)  = X_next(1:num_masses, 1);
            obj.vrx(~obj.slide)      = X_next(num_masses+1:2 * num_masses, 1);
            obj.delta_y(~obj.slide)  = X_next(2 * num_masses+1: 3 * num_masses, 1);
            obj.vry(~obj.slide)      = X_next(3 * num_masses+1:4 * num_masses, 1); 

           % % Update velocity history
           %  obj.prev3_vrx(~obj.slide)   = obj.prev2_vrx(~obj.slide);
           %  obj.prev3_vry(~obj.slide)   = obj.prev2_vry(~obj.slide);
           %  obj.prev2_vrx(~obj.slide)   = obj.prev1_vrx(~obj.slide);
           %  obj.prev2_vry(~obj.slide)   = obj.prev1_vry(~obj.slide);
           %  obj.prev1_vrx(~obj.slide)   = obj.vrx(~obj.slide);
           %  obj.prev1_vry(~obj.slide)   = obj.vry(~obj.slide);
           % 
           %  % Calculate New Velocities
           %  obj.vrx(~obj.slide) = omega .* re + omega_z .* (obj.y(~obj.slide) + obj.delta_y(~obj.slide)) - v0 .* cos(alpha);
           %  obj.vry(~obj.slide) = -omega_z .* (obj.x(~obj.slide) + obj.delta_x(~obj.slide)) - v0 .* sin(alpha);
           % 
           %  % 2nd order backward difference for acceleration estimation
           %  obj.dvrx(~obj.slide) = (3 .* obj.vrx(~obj.slide) - 4 .* obj.prev1_vrx(~obj.slide) + obj.prev2_vrx(~obj.slide)) ./ (2 * dt);
           %  obj.dvry(~obj.slide) = (3 .* obj.vry(~obj.slide) - 4 .* obj.prev1_vry(~obj.slide) + obj.prev2_vry(~obj.slide)) ./ (2 * dt);
           % 
            % % Update displacements for sticking using verlet integration
            % obj.delta_x(~obj.slide) = 2 * obj.delta_x(~obj.slide) - obj.prev1_delta_x(~obj.slide) + obj.dvrx(~obj.slide) * dt^2;
            % obj.delta_y(~obj.slide) = 2 * obj.delta_y(~obj.slide) - obj.prev1_delta_y(~obj.slide) + obj.dvrx(~obj.slide) * dt^2;

            % Calculate shear stres
            obj.tauX(~obj.slide) = obj.kx .* obj.delta_x(~obj.slide) + obj.cx .* obj.vrx(~obj.slide) + obj.m .* obj.dvrx(~obj.slide);
            obj.tauY(~obj.slide) = obj.ky .* obj.delta_y(~obj.slide) + obj.cy .* obj.vry(~obj.slide) + obj.m .* obj.dvry(~obj.slide);
            
            % % Update displacement history
            % obj.prev1_delta_x = obj.delta_x;
            % obj.prev1_delta_y = obj.delta_y;
        end

        function [obj] = update_brush(obj, pressVal, omega, omega_z, re, v0, alpha, dt, t)
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
                t           (1,1)   single
            end
            
            % Validate inputs
            validateattributes(pressVal, {'single'}, {'real'})
            validateattributes(omega, {'single'}, {'scalar'});
            validateattributes(omega_z, {'single'}, {'scalar'});
            validateattributes(re, {'single'}, {'scalar'});
            validateattributes(v0, {'single'}, {'scalar'});
            validateattributes(alpha, {'single'}, {'scalar'});
            validateattributes(dt, {'single'}, {'scalar', 'positive'});
            validateattributes(t, {'single'}, {'scalar', 'nonnegative'});

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
                        slide_obj = slide_obj.solve_brush_slide_ode(omega_vec(obj.slide), omega_z, re, v0_vec(obj.slide), alpha, dt, t);
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
                    obj = obj.solve_brush_slide_ode(omega_vec, omega_z, re, v0_vec, alpha, dt, t);
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

        function obj = calculateStresses(obj,~, X,forcing_args)
            %CALCULATESTRESSES a method used in the time stepping of the
            %dynamical system
            % INPUTS
            % obj:      Brush object containing the previous properties
            % X:        
            % omega:    
            % re:       
            % omega_z:  
            % v0:       
            % alpha:    
            %
            % OUTPUTS
            % obj:      Brush object containing the new properties
            
            % Extract additional parameters
            omega   = forcing_args.omega;
            re      = forcing_args.re;
            omega_z = forcing_args.omega_z;
            v0      = forcing_args.v0;
            alpha   = forcing_args.alpha;
            
            % Extract states from state vector
            num_masses = length(obj.vrx(obj.slide));
            DeltaX = X(1:num_masses);
            DeltaY = X(num_masses + 1:2 * num_masses);
            VX     = X(2 * num_masses + 1:3 * num_masses);
            VY     = X(3 * num_masses + 1:4 * num_masses);
            
            % -------------------------------------------------------------
            %               Main Calculations start here
            % -------------------------------------------------------------
            
            obj.vrx(obj.slide) = omega .* re + omega_z .* (obj.y(obj.slide) + DeltaY) - v0 .* cos(alpha);
            obj.vry(obj.slide) = -omega_z .* (obj.x(obj.slide) + DeltaX) - v0 .* sin(alpha);
            
            % Calculate sliding velocity
            obj.vs_y(obj.slide) = ( VY - obj.vry(obj.slide) );
            obj.vs_x(obj.slide) = ( VX - obj.vrx(obj.slide) );
            obj.vs(obj.slide) = ( hypot(obj.vs_x(obj.slide), obj.vs_y(obj.slide)) ) .* sign(obj.vs_x(obj.slide));

            obj.theta_2(obj.slide) = real( atan2(obj.vs_y(obj.slide), obj.vs_x(obj.slide)) ) - pi; % Minus pi to enforce colinearity of stress and velocity vectors

            % Calculate friction coefficient
            obj.mu(obj.slide) = obj.compute_friction_coefficient(obj.press(obj.slide), obj.p_0, obj.p_ref,...
                                                                 obj.q, obj.mu_0, obj.mu_m, obj.h, obj.vs(obj.slide), obj.v_m);

            % Calculate stresses
            obj.tauX(obj.slide) = obj.mu(obj.slide) .* obj.press(obj.slide) .* cos(obj.theta_2(obj.slide));
            obj.tauY(obj.slide) = obj.mu(obj.slide) .* obj.press(obj.slide) .* sin(obj.theta_2(obj.slide));
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
            slide = (tau > mu_0 .* press);
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
                press   (:,1) single
                p_0     (1,1)   single
                p_ref   (1,1)   single
                q       (1,1)   single
                mu_0    (1,1)   single
                mu_m    (1,1)   single
                h       (1,1)   single
                vs      (:,1) single
                v_m     (1,1)   single
            end
            
            % Calculate pressure-dependent term
            masked_p = max(press, p_0);
            sat_p = real( (masked_p ./ p_ref) .^(-q) );

            % Ensure velocity is valid for logarithm
            vm_safe = max(abs(v_m), eps);
            vs_safe = max(abs(vs), eps);
            
            % Calculate velocity-dependent term
            logterm = log10(vs_safe ./ vm_safe);
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

            