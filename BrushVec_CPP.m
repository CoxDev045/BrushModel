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
        kx          (1,1) single = 10;            % Base x-stiffness
        ky          (1,1) single = 10;%0.37;      % Base y-stiffness
        cx          (1,1) single = 0.1;%1.78e-7;   % x-damping coefficient
        cy          (1,1) single = 0.1;%1.40e-4;   % y-damping coefficient
        m_inv       (1,1) single = 10;%1.3089005e+09;  % Inverse of Mass property
        m           (1,1) single = 0.1;%7.64e-10;  % Mass

        % Friction Model Properties (Static)
        mu_0        (1,1) single = 0.02;      % Static friction coefficient
        mu_m        (1,1) single = 1.15;      % Maximum friction coefficient
        h           (1,1) single = 0.4;%0.23;      % Friction model parameter
        p_0         (1,1) single = 0.02;      % Minimum pressure threshold
        p_ref       (1,1) single = 0.39;      % Reference pressure
        q           (1,1) single = 0.28;      % Pressure exponent
        v_m         (1,1) single = 5;%23.47;     % Reference velocity

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

            % Initialize properties
            obj.x = xVal;
            obj.minX = min(xVal);
            obj.maxX = max(xVal);
            obj.y = yVal;
            obj.press = press;                       % Vertical Pressure

            initConds = zeros(size(xVal, 1), 22);
            initConds(:,1) = obj.mu_0;
            initConds_bool = false(size(xVal, 1), 2);

            obj = obj.initialConditions(initConds, initConds_bool);

        end

        function [obj] = initialConditions(obj, initConds, initConds_bool)
            
            obj.slide = initConds_bool(:, 1);
            obj.passed = initConds_bool(:, 2);

            obj.mu = initConds(:, 1);
            
            obj.vs = initConds(:, 2);
            
            obj.delta_x = initConds(:, 3);
            obj.delta_y = initConds(:, 4);
            
            obj.tauX = initConds(:, 5);
            obj.tauY = initConds(:, 6);

            obj.vrx = initConds(:, 7);
            obj.vry = initConds(:, 8);

            obj.prev3_vrx = initConds(:, 9);          % Relative x velocity, 3 steps back
            obj.prev3_vry = initConds(:, 10);         % Relative y velocity, 3 steps back
            obj.prev2_vrx = initConds(:, 11);         % Relative x velocity, 2 steps back
            obj.prev2_vry = initConds(:, 12);         % Relative y velocity, 2 steps back
            obj.prev1_vrx = initConds(:, 13);         % Relative x velocity, 1 step back
            obj.prev1_vry = initConds(:, 14);         % Relative y velocity, 1 step back
            obj.theta_2 = initConds(:, 15);           % Angle between sliding velocity and horizontal
            obj.theta_1 = initConds(:, 16);           % Angle between shear stress and horizontal

            % Integration method state variables
            obj.ax = initConds(:, 17);         % x-acceleration
            obj.ay = initConds(:, 18);         % y-acceleration
            obj.ax_new = initConds(:, 19);     % New x-acceleration for Verlet
            obj.ay_new = initConds(:, 20);     % New y-acceleration for Verlet
            obj.vx = initConds(:, 21);         % x-velocity
            obj.vy = initConds(:, 22);         % y-velocity
        end

        function [finalProps, finalSlide] = saveFinalProperties(obj)
            
            finalProps = zeros(400, 22);
            finalSlide = false(400,1);

            finalProps(:, 1) = obj.mu;
            
            finalProps(:, 2) = obj.vs;

            finalSlide = obj.slide;
            
            finalProps(:, 3) = obj.delta_x;
            finalProps(:, 4) = obj.delta_y;
           
            finalProps(:, 5) = obj.tauX;
            finalProps(:, 6) = obj.tauY;

            finalProps(:, 7) = obj.vrx;
            finalProps(:, 8) = obj.vry;

            finalProps(:, 9) = obj.prev3_vrx;          % Relative x velocity, 3 steps back
            finalProps(:, 10) = obj.prev3_vry;         % Relative y velocity, 3 steps back
            finalProps(:, 11) = obj.prev2_vrx;         % Relative x velocity, 2 steps back
            finalProps(:, 12) = obj.prev2_vry;         % Relative y velocity, 2 steps back
            finalProps(:, 13) = obj.prev1_vrx;         % Relative x velocity, 1 step back
            finalProps(:, 14) = obj.prev1_vry;         % Relative y velocity, 1 step back
            finalProps(:, 15) = obj.theta_2;           % Angle between sliding velocity and horizontal
            finalProps(:, 16) = obj.theta_1;           % Angle between shear stress and horizontal

            % Integration method state variables
            finalProps(:, 17) = obj.ax;         % x-acceleration
            finalProps(:, 18) = obj.ay;         % y-acceleration
            finalProps(:, 19) = obj.ax_new;     % New x-acceleration for Verlet
            finalProps(:, 20) = obj.ay_new;     % New y-acceleration for Verlet
            finalProps(:, 21) = obj.vx;         % x-velocity
            finalProps(:, 22) = obj.vy;         % y-velocity
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

            obj.vrx(obj.slide) = omega .* re + omega_z .* (obj.y(obj.slide) + obj.delta_y(obj.slide)) - v0 .* cos(alpha);
            obj.vry(obj.slide) = -omega_z .* (obj.x(obj.slide) + obj.delta_x(obj.slide)) - v0 .* sin(alpha);
            
            % Calculate sliding velocity
            obj.vs_y(obj.slide) = real( obj.vy(obj.slide) - obj.vry(obj.slide) );
            obj.vs_x(obj.slide) = real( obj.vx(obj.slide) - obj.vrx(obj.slide) );
            obj.vs(obj.slide) = real( hypot(obj.vs_x(obj.slide), obj.vs_y(obj.slide)) );

            obj.theta_1(obj.slide) = real( atan2(obj.vs_y(obj.slide), obj.vs_x(obj.slide)) ); 
            obj.theta_2(obj.slide) = obj.theta_1(obj.slide) - pi;

            % % cos_theta_1 = cos(obj.theta_1(obj.slide));
            % % mask = abs(cos_theta_1) > 1e-8;
            % % obj.vs(obj.slide) = mask .* (obj.vx(obj.slide) - obj.vrx(obj.slide)) ./ cos_theta_1;

            % Calculate friction coefficient
            obj.mu = obj.compute_friction_coefficient(obj.press, obj.p_0, obj.p_ref, obj.q, obj.mu_0, obj.mu_m, obj.h, obj.vs, obj.v_m);

            % Calculate stresses
            obj.tauX(obj.slide) = obj.mu(obj.slide) .* obj.press(obj.slide) .* cos(obj.theta_2(obj.slide));
            obj.tauY(obj.slide) = obj.mu(obj.slide) .* obj.press(obj.slide) .* sin(obj.theta_2(obj.slide));

            % Apply selected integration method
            obj = obj.integrate_verlet(dt);
            
        end
        
        function [obj] = integrate_verlet(obj, dt)
            % Verlet integration method for brush dynamics
            %
            % Parameters:
            %   dt - Time step
            
            % Update Velocity History
            obj.prev1_vx = obj.vx;
            obj.prev1_vy = obj.vy;
            
            % First half-step velocity update
            obj.vx(obj.slide) = obj.prev1_vx(obj.slide) + 0.5 .* obj.ax(obj.slide) .* dt;
            obj.vy(obj.slide) = obj.prev1_vy(obj.slide) + 0.5 .* obj.ay(obj.slide) .* dt;
            
            % Full position update
            obj.delta_x(obj.slide) = obj.delta_x(obj.slide) + obj.vx(obj.slide) .* dt + 0.5 .* obj.ax(obj.slide) .* dt^2;
            obj.delta_y(obj.slide) = obj.delta_y(obj.slide) + obj.vy(obj.slide) .* dt + 0.5 .* obj.ay(obj.slide) .* dt^2;
            
            % Calculate new accelerations
            obj.ax_new(obj.slide) = (obj.tauX(obj.slide) - obj.kx .* obj.delta_x(obj.slide) - obj.cx .* obj.vx(obj.slide)) * obj.m_inv;
            obj.ay_new(obj.slide) = (obj.tauY(obj.slide) - obj.ky .* obj.delta_y(obj.slide) - obj.cy .* obj.vy(obj.slide)) * obj.m_inv;
            
            % Second half-step velocity update with new accelerations
            obj.vx(obj.slide) = obj.vx(obj.slide) + 0.5 .* (obj.ax_new(obj.slide) + obj.ax(obj.slide)) .* dt;
            obj.vy(obj.slide) = obj.vy(obj.slide) + 0.5 .* (obj.ay_new(obj.slide) + obj.ay(obj.slide)) .* dt;
            
            % Update accelerations for next step
            obj.ax(obj.slide) = obj.ax_new(obj.slide);
            obj.ay(obj.slide) = obj.ay_new(obj.slide);
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
            
            % % p0_vec = p_0 * ones(size(press));
            % % p_ref_vec = p_ref * ones(size(press));
            % q_vec = q * ones(size(press));

            % Calculate pressure-dependent term
            % % mask_p = press > p_0;
            % % sat_p = zeros(size(press));
            % % sat_p(mask_p) = real( (press(mask_p) ./ p_ref_vec(mask_p)).^(-q) );
            % % sat_p(~mask_p) = real( (p0_vec(~mask_p) ./ p_ref_vec(~mask_p)).^(-q) );

            % % masked_p = max(press, p_0);
            sat_p = real( (press ./ p_ref) .^(-q) ); % .* mask_p + ((p_0 / p_ref) .^ (-q)) .* ~mask_p;

            % % clear p0_vec p_ref_vec q_vec mask_p;

            % Ensure velocity is valid for logarithm
            vm_safe = max(abs(v_m), eps) .* sign(v_m);
            vs_safe = max(abs(vs), eps) .* sign(vs);
            
            % Calculate velocity-dependent term
            logterm = log10(abs(vs_safe ./ vm_safe));
            scaleTerm = (mu_m - mu_0);
            expTerm = exp(-h.^2 .* logterm.^2);
            mu = sat_p .* (mu_0 +  scaleTerm .* expTerm);

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

            