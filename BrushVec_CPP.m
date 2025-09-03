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
        mu_0        (1,1) single = 0.02;                        % Static friction coefficient
        mu_m        (1,1) single = 1.20;                        % Maximum friction coefficient
        h           (1,1) single = 1.5;%0.4;                    % Friction model parameter
        p_0         (1,1) single = 0.02;                        % Minimum pressure threshold
        p_ref       (1,1) single = 0.247645523118594;%0.39;     % Reference pressure
        q           (1,1) single = 0.390845735345209;%0.28;     % Pressure exponent
        v_m         (1,1) single = 5;%30.206780784527050;%23.47;% Reference velocity

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
        mu              (400,1) single  = 0.02;       % Friction coefficient
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
            forcing_args = struct('SlideInd', slideInd ,...
                                  'x', obj.x(slideInd),... 
                                  'y', obj.y(slideInd),... 
                                  'omega', omega,... 
                                  're', re,... 
                                  'omega_z', omega_z,... 
                                  'v0', v0,... 
                                  'alpha', alpha,... 
                                  'press', obj.press(slideInd),... 
                                  'p_0', obj.p_0,... 
                                  'p_ref', obj.p_ref,... 
                                  'q', obj.q,... 
                                  'mu_0', obj.mu_0,... 
                                  'mu_m', obj.mu_m,... 
                                  'h', obj.h,... 
                                  'v_m', obj.v_m);

            brush_dynamics = @(t, X) BrushVec_CPP.slidingDynamics(t,X, forcing_args, obj.m, obj.kx, obj.ky, obj.cx, obj.cy);
            args = cell(1);
            X_next = trbdf2_step(brush_dynamics, dt, t, X_vec, args);
            
            %
            num_masses = length(X_next) / 4;
            % Extract into object
            obj.delta_x(slideInd) = X_next(1:num_masses, 1);
            obj.vx(slideInd)      = X_next(num_masses+1:2 * num_masses, 1);
            obj.delta_y(slideInd) = X_next(2 * num_masses+1: 3 * num_masses, 1);
            obj.vy(slideInd)      = X_next(3 * num_masses+1:4 * num_masses, 1); 
           
        end

        function [obj] = solve_stick_ode(obj, omega, omega_z, re, v0, alpha, dt)
            %SOLVE_STICK_ODE is a method that steps the simulation forward
            %by one time step. This method is only called when there are
            %brushes that are not sliding. This means there is only
            %relative acceleration between road element and brush due to
            %the relative velocity
            % (d^δ/dt^2 = d^2v_r/dt
            % => 𝜏 = k * δ + c * dv_r/dt + m * d^2v_r/dt)
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

            % Update displacement history
            obj.prev1_delta_x = obj.delta_x;
            obj.prev1_delta_y = obj.delta_y;
        
            % Calculate shear stres
            obj.tauX(~obj.slide) = obj.kx .* obj.delta_x(~obj.slide) + obj.cx .* obj.vrx(~obj.slide) + obj.m .* obj.dvrx(~obj.slide);
            obj.tauY(~obj.slide) = obj.ky .* obj.delta_y(~obj.slide) + obj.cy .* obj.vry(~obj.slide) + obj.m .* obj.dvry(~obj.slide);
        
            % Calculate New Velocities
            obj.vrx(~obj.slide) = omega .* re + omega_z .* (obj.y(~obj.slide) + obj.delta_y(~obj.slide)) - v0 .* cos(alpha);
            obj.vry(~obj.slide) = -omega_z .* (obj.x(~obj.slide) + obj.delta_x(~obj.slide)) - v0 .* sin(alpha);

            % 2nd order backward difference for acceleration estimation
            obj.dvrx(~obj.slide) = (3 .* obj.vrx(~obj.slide) - 4 .* obj.prev1_vrx(~obj.slide) + obj.prev2_vrx(~obj.slide)) ./ (2 * dt);
            obj.dvry(~obj.slide) = (3 .* obj.vry(~obj.slide) - 4 .* obj.prev1_vry(~obj.slide) + obj.prev2_vry(~obj.slide)) ./ (2 * dt);

            % Update displacements for sticking
            obj.delta_x(~obj.slide) = 2 * obj.delta_x(~obj.slide) - obj.prev1_delta_x(~obj.slide) + obj.dvrx(~obj.slide) * dt^2;
            obj.delta_y(~obj.slide) = 2 * obj.delta_y(~obj.slide) - obj.prev1_delta_y(~obj.slide) + obj.dvrx(~obj.slide) * dt^2;
            
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
            validateattributes(omega, {'single'}, {'scalar'});
            validateattributes(omega_z, {'single'}, {'scalar'});
            validateattributes(re, {'single'}, {'scalar'});
            validateattributes(v0, {'single'}, {'scalar'});
            validateattributes(alpha, {'single'}, {'scalar'});
            validateattributes(dt, {'single'}, {'scalar', 'positive'});
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
        function [tauX, tauY] = calculateStresses(t, X, args)
            %CALCULATESTRESSES a method used in the time stepping of the
            %dynamical system
            % INPUTS
            % X:        Current state vector containing the following:
            %           [delta_x(:), delta_y(:), vx(:), vy(:)];
            %
            % args:     Struct containing the following:     
            %           SlideInd:   Vector containing the indices where sliding occurs
            %           x:          X-Coordinates of brushes
            %           y:          Y-coordinates of brushes
            %           omega:      Current angular velocity
            %           re:         Effective rolling radius
            %           omega_z:    Turning velocity of wheel
            %           v0:         Velocity of road element
            %           alpha:      Slip angle
            %           press:      Current vertical pressure value
            %           p_0:        Minimum pressure value allowed
            %           p_ref:      Reference pressure value
            %           q:          Pressure scaling value
            %           mu_0:       Minimum friction coefficient
            %           mu_m:       Maximum friction coefficient
            %           h:          Friction scaling parameter
            %           v_m:        Speed where max friction coefficient occurs
            %           
            % OUTPUTS
            % tauX:     Vector containing the longitudinal stresses based on the current
            %           state estimate
            % tauY:     Vector containing the lateral stresses based on the current
            %           state estimate
            
            % The total number of masses is the length of kx
            num_masses = length(X) / 4; % 4 states present so will always be a factor of 4
        
            % Unpack the state vector X (a single column vector)
            delta_x = X(1:num_masses);
            delta_y = X(num_masses+1 : 2*num_masses);
            vx = X(2*num_masses+1 : 3*num_masses);
            vy = X(3*num_masses+1 : 4*num_masses);
        
            % Initialise intermediate matrices
            % vrx = zeros(size(vx));
            % vry = zeros(size(vy));
            % vs_x = zeros(size(vx));
            % vs_y = zeros(size(vy));
            % vs = zeros(size(vx));
            % theta_2 = zeros(size(vx));
            % tauX = zeros(size(vx));
            % tauY = zeros(size(vx));
            
            % Extract additional parameters
            % slideInd = args.SlideInd;
            x = args.x; y = args.y;
            omega = args.omega; re = args.re;
            omega_z = args.omega_z;
            v0 = args.v0; alpha = args.alpha;
            press = args.press;
            p_0 = args.p_0; p_ref = args.p_ref;
            q = args.q;
            mu_0 = args.mu_0; mu_m = args.mu_m;
            h = args.h;
            v_m = args.v_m;
            
            % -------------------------------------------------------------
            %               Main Calculations start here
            % -------------------------------------------------------------
            vrx = omega .* re + omega_z .* (y + delta_y) - v0 .* cos(alpha);
            vry = -omega_z .* (x + delta_x) - v0 .* sin(alpha);
            
            % Calculate sliding velocity
            vs_y = ( vy - vry );
            vs_x = ( vx - vrx );
            vs = ( hypot(vs_x, vs_y) ) .* sign(vs_x);

            theta_2 = real( atan2(vs_y, vs_x) ) - pi; % Minus pi to enforce colinearity of stress and velocity vectors

            % Calculate friction coefficient
            mu = BrushVec_CPP.compute_friction_coefficient(press, p_0, p_ref, q, mu_0, mu_m, h, vs, v_m);

            % Calculate stresses
            tauX = mu .* press .* cos(theta_2);
            tauY = mu .* press .* sin(theta_2);

        end

        function dX = slidingDynamics(t,X, forcing_args, m, kx, ky, cx, cy)
            % INPUTS
            % X:                Current state vector containing the following:
            %                   [delta_x(:), delta_y(:), vx(:), vy(:)];
            % forcing_func:     Function handle of forcing function
            % forcing_args:     Struct containing the parameters for the force calculation
            %       SlideInd:   Vector containing the indices where sliding occurs
            %       x:          X-Coordinates of brushes
            %       y:          Y-coordinates of brushes
            %       omega:      Current angular velocity
            %       re:         Effective rolling radius
            %       omega_z:    Turning velocity of wheel
            %       v0:         Velocity of road element
            %       alpha:      Slip angle
            %       press:      Current vertical pressure value
            %       p_0:        Minimum pressure value allowed
            %       p_ref:      Reference pressure value
            %       q:          Pressure scaling value
            %       mu_0:       Minimum friction coefficient
            %       mu_m:       Maximum friction coefficient
            %       h:          Friction scaling parameter
            %       v_m:        Speed where max friction coefficient occurs
            %
            % m:                System inertia parameter
            % kx / ky:          System stiffness coefficients (in x and y directions)
            % cx / cy:          System damping coefficients (in x and y directions)
            % OUTPUTS
            % dx:       Vector containing the derivatives of the next state
            %           [ddelta_x;
            %            ddelta_y;
            %            dvx;
            %            dvy];
            %           
            
            SlideInd = forcing_args.SlideInd;
            % The total number of masses is the length of kx
            num_masses = length(X) / 4; % 4 states present so will always be a factor of 4
        
            % Unpack the state vector X (a single column vector)
            delta_x = X(1:num_masses);
            delta_y = X(num_masses+1 : 2*num_masses);
            vx =      X(2*num_masses+1 : 3*num_masses);
            vy =      X(3*num_masses+1 : 4*num_masses);
        
            % Initialize derivative vectors
            d_delta_x = zeros(num_masses, 1);
            d_vx = zeros(num_masses, 1);
            d_delta_y = zeros(num_masses, 1);
            d_vy = zeros(num_masses, 1);
            
            % Get forcing at time t
            [Fx, Fy] = BrushVec_CPP.calculateStresses(t, X, forcing_args);
        
            % Apply dynamics only where sliding occurs
            d_delta_x = vx;
            d_vx      = (Fx - kx .* delta_x - cx .* vx) ./ m;
        
            d_delta_y = vy;
            d_vy      = (Fy - ky .* delta_y - cy .* vy) ./ m;
        
            % Reassemble the output into a single column vector
            dX = [d_delta_x; d_delta_y; d_vx; d_vy];           
        end

        function dX = stickingDynamics(t,X, SlideInd, m, kx, ky, cx, cy, Forcing)
            % INPUTS
            % X:        Current state vector containing the following:
            %           [delta_x(:), delta_y(:), vx(:), vy(:)];
            % SlideInd: Vector containing the indices where sliding occurs
            % m:        
            % kx / ky:  System stiffness coefficients (in x and y directions)
            % cx / cy:  System damping coefficients (in x and y directions)
            % F:        Vector containing the stresses which act as a forcing
            %           function: [tauX(:), tauY(:)];
            % OUTPUTS
            % dx:       Vector containing the derivatives of the next state
            %           [dX(:);
            %            dY(:);
            %           
            
            % The total number of masses is the length of kx
            num_masses = length(kx);
        
            % Unpack the state vector X (a single column vector)
            delta_x = X(1:num_masses);
            delta_y = X(num_masses+1 : 2*num_masses);
            vx = X(2*num_masses+1 : 3*num_masses);
            vy = X(3*num_masses+1 : 4*num_masses);
        
            % Initialize derivative vectors
            d_delta_x = zeros(num_masses, 1);
            d_vx = zeros(num_masses, 1);
            d_delta_y = zeros(num_masses, 1);
            d_vy = zeros(num_masses, 1);
            
            % Get forcing at time t
            [Fx, Fy] = Forcing(t, X);
        
            % Apply dynamics only where sliding occurs
            d_delta_x(SlideInd) = vx(SlideInd);
            d_vx(SlideInd) = (Fx(SlideInd) - kx(SlideInd) .* delta_x(SlideInd) - cx(SlideInd) .* vx(SlideInd)) ./ m(SlideInd);
        
            d_delta_y(SlideInd) = vy(SlideInd);
            d_vy(SlideInd) = (Fy(SlideInd) - ky(SlideInd) .* delta_y(SlideInd) - cy(SlideInd) .* vy(SlideInd)) ./ m(SlideInd);
        
            % Reassemble the output into a single column vector
            dX = [d_delta_x; d_delta_y; d_vx; d_vy];           
        end

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

            