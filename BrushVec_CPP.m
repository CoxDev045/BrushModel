classdef BrushVec_CPP < handle%#codegen -args
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
        kx          (1,1) single = 2.5e-4;            % Base x-stiffness
        ky          (1,1) single = 2.5e-4;%0.37;      % Base y-stiffness
        kz          (:,1) single = 0; % Area normalised stiffness (N/m) /m^2
        kz_base     (1,1) single = 200e3% Area normalised base stiffness ( N/(m * sqrt(Pa)) ) /m^2
        cx          (1,1) single = 1.5e-7;%1.78e-7;   % x-damping coefficient
        cy          (1,1) single = 1.5e-4;%1.40e-4;   % y-damping coefficient
        cz          (1,1) single = 250% Area normalised damping   (Ns/m)/m^2
        m           (1,1) single = 1e-10;%0.00106768549514851;%7.64e-10;  % Mass

        % Friction Model Properties (Static)
        mu_0        (1,1) double = 0.008918;        % Static friction coefficient
        mu_m        (1,1) double = 0.941811;        % Maximum friction coefficient
        h           (1,1) double = 0.638421;        % Friction model shape parameter
        p_0         (1,1) double = 0.02;            % Minimum pressure threshold
        p_ref       (1,1) double = 0.119023;        % Reference pressure
        q           (1,1) double = 0;%0.019526;       % Pressure exponent
        v_m         (1,1) double = 4.248101;        % Reference velocity

        % Dynamic properties (Changing throughout simulation)
        numBrushes      (1, 1) uint16   = 400;       % Number of brushes in model
        minX            (1, 1)  single  = -1;       % Minimum x-value that x can be in when brush is in contact patch
        maxX            (1, 1)  single  = 1;       % Maximum x_value that x can be when brush is in contact patch
        x               (:,1) single  = 0;       % x-coordinate of the brush
        y               (:,1) single  = 0;       % y-coordinate of the brush
        delta_x         (:,1) single  = 0;       % x-displacement of the brush
        delta_y         (:,1) single  = 0;       % y-displacement of the brush
        tauX            (:,1) single  = 0;       % x shear force
        tauY            (:,1) single  = 0;       % y shear force
        mu              (:,1) single  = 0.02;       % Friction coefficient
        press           (:,1) single  = 0;       % Current Pressure
        initPress       (:,1) single  = 0;       % Initial Pressure
        vrx             (:,1) single  = 0;       % Relative x velocity
        vry             (:,1) single  = 0;       % Relative y velocity
        vs              (:,1) single  = 0;       % Sliding velocity magnitude
        vs_x            (:,1) single  = 0;       % X Sliding Velocity
        vs_y            (:,1) single  = 0;       % Y Sliding Velocity
        theta_2         (:,1) single  = 0;       % Angle between sliding velocity and horizontal
        slide           (:,1) logical = false;       % Sliding state (true/false)
        hasPress        (:,1) logical = false;
        passed          (:,1) logical = false;       % To check whether the brush has moved past contact length
        vx              (:,1) single = 0;        % X Deformation Velocity
        vy              (:,1) single = 0;        % Y Deformation Velocity

        X_next          (:, 1) single = 0;

    end
    
    methods
        % Constructor
        function obj = BrushVec_CPP(xVal, yVal, press, numBrushes)
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

            % % Find these parameters seperately from data
            % obj.h = 0.75;%2.0;     % Stribeck Exponent
            % obj.q = 0.390845735345209;     % Pressure Saturation Exponent
            % obj.v_m = 2.5;%30.206780784527050;   % Stribeck Velocity

            
            % Initialize properties from inputs
            obj.x = xVal;
            obj.minX = min(xVal);
            obj.maxX = max(xVal);
            obj.y = yVal;
            obj.press = press;         % Current Vertical Pressure
            obj.initPress = press;     % Initial vertical Pressure

            obj.kz_base = obj.kz_base / single(obj.numBrushes);
            obj.cz = obj.cz / single(obj.numBrushes);
            

            obj.kz = obj.kz_base .* (max(obj.press, 0).^(0.5) );
            % Initialise other properties that are not zero
            obj.mu = obj.mu_0;

            % Initialise other properties that are zero
            obj.delta_x = zeros(size(press));
            obj.delta_y = zeros(size(press));
            obj.vrx = zeros(size(press));
            obj.vry = zeros(size(press));
            obj.vs_x = zeros(size(press));
            obj.vs_y = zeros(size(press));
            obj.vs = zeros(size(press));            
            obj.vx = zeros(size(press));
            obj.vy = zeros(size(press));
            obj.tauX = zeros(size(press));
            obj.tauY = zeros(size(press));
            obj.theta_2 = zeros(size(press));
            obj.mu = zeros(size(press));
            
            obj.X_next = zeros(4 * obj.numBrushes, 1);
            
            

        end

        function obj = update_brush(obj, pressVal, omega, omega_z, re, v0, alpha, dt, t, z_disp, z_dot)
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
                pressVal    (:,1) single
                omega       
                omega_z     (1,1)   single
                re          (1,1)   single
                v0          
                alpha       (1,1)   single
                dt          (1,1)   single
                t           (1,1)   single
                z_disp      (:,1)   single
                z_dot       (:,1)   single
            end
            
            % Validate inputs
            validateattributes(pressVal, {'single'}, {'real'})
            % validateattributes(omega, {'single'}, {'scalar'});
            validateattributes(omega_z, {'single'}, {'scalar'});
            validateattributes(re, {'single'}, {'scalar'});
            % validateattributes(v0, {'single'}, {'scalar'});
            validateattributes(alpha, {'single'}, {'scalar'});
            validateattributes(dt, {'single'}, {'scalar', 'positive'});
            validateattributes(t, {'single'}, {'scalar', 'nonnegative'});

            % % % Assume uniform angular velocity and translational velocity of
            % % % road element is constant across contact patch
            % % omega_vec = repmat(omega, size(pressVal));
            % % v0_vec = repmat(v0, size(pressVal));

            %%%%%%%%%%%%%% Update Pressure Dependent Properties %%%%%%%%%%%%%%%%
            obj.press(1:obj.numBrushes) = pressVal(1:obj.numBrushes);
            % Determine whether elements have vertical pressure applied
            obj.hasPress = obj.has_pressure(obj.press, obj.p_0);
            % Update Pressure distribution
            obj = obj.calculateVertical(z_disp, z_dot);
            %%%%%%%%%%%%% Remove Pressure Dependent stiffness for now %%%%%
            % % obj.ky = obj.ky_0 + obj.ky_l .* obj.press;
            % % obj.kx = obj.phi .* obj.ky;
            
            % Determine if elements are sliding or sticking
            obj.slide = obj.is_sliding(obj.tauX, obj.tauY, obj.mu_0, obj.press);

            % Solve dynamics for entire system
            obj.solve_dynamics_ode(omega, omega_z, re, v0, alpha, dt, t);
            
            % % Update dynamics based on state
            % if any(obj.slide)
            %     % Handle partial sliding
            %     if ~all(obj.slide)
            %         % Split object to handle slide and stick separately
            % 
            %         if ~isempty(omega_vec(obj.slide))
            %             slide_obj = obj;
            %             % Process sliding elements
            %             slide_obj = slide_obj.solve_brush_slide_ode(omega_vec(obj.slide), omega_z, re, v0_vec(obj.slide), alpha, dt, t);
            %             % Combine results back
            %             obj.delta_x(obj.slide) = slide_obj.delta_x(obj.slide);
            %             obj.delta_y(obj.slide) = slide_obj.delta_y(obj.slide);
            %             obj.vx(obj.slide) = slide_obj.vx(obj.slide);
            %             obj.vy(obj.slide) = slide_obj.vy(obj.slide);
            % 
            %             obj.vrx(obj.slide) = slide_obj.vrx(obj.slide);
            %             obj.vry(obj.slide) = slide_obj.vry(obj.slide);
            % 
            %             obj.tauX(obj.slide) = slide_obj.tauX(obj.slide);
            %             obj.tauY(obj.slide) = slide_obj.tauY(obj.slide);
            %         end
            % 
            %         if ~isempty(omega_vec(~obj.slide))
            %             stick_obj = obj;
            %             % Process sticking elements
            %             stick_obj = stick_obj.solve_stick_ode(omega_vec(~obj.slide), omega_z, re, v0_vec(~obj.slide), alpha, dt);
            %             % Combine results back
            %             obj.delta_x(~obj.slide) = stick_obj.delta_x(~obj.slide);
            %             obj.delta_y(~obj.slide) = stick_obj.delta_y(~obj.slide);
            % 
            %             obj.vrx(~obj.slide) = stick_obj.vrx(~obj.slide);
            %             obj.vry(~obj.slide) = stick_obj.vry(~obj.slide);
            % 
            %             obj.tauX(~obj.slide) = stick_obj.tauX(~obj.slide);
            %             obj.tauY(~obj.slide) = stick_obj.tauY(~obj.slide);
            %         end
            %     else
            %         % All elements are sliding
            %         obj = obj.solve_brush_slide_ode(omega_vec, omega_z, re, v0_vec, alpha, dt, t);
            %     end
            % else
            %     % All elements are sticking
            %     obj = obj.solve_stick_ode(omega_vec, omega_z, re, v0_vec, alpha, dt);
            % end
            
            obj.updateXvalues(omega(t), re, dt);
            
        end

        function obj = solve_dynamics_ode(obj, omega, omega_z, re, v0, alpha, dt, t)
            arguments
                obj
                omega   
                omega_z (1,1) single
                re      (1,1) single
                v0      
                alpha   (1,1) single
                dt      (1,1) single
                t       (1,1) single
            end
            % - get current state vector of X = [delta_x; delta_y;vx;vy];
            X_vec = [obj.delta_x; obj.delta_y; obj.vx; obj.vy];

	        % - pass state vector to integrator: X_next = integrateDynamics(t, X, dt) where some integration scheme is applied to step X to X_next via dX
            obj = obj.integrateDynamics(dt, t, X_vec, 'rk4', omega, omega_z, re, v0, alpha);
	
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
                obj.tauX(hasPassed) = obj.mu_0 .* obj.press(hasPassed) .* cos(-pi);
                obj.tauY(hasPassed) = obj.mu_0 .* obj.press(hasPassed) .* sin(-pi);
                %%% Reset friction coefficient
                obj.mu(hasPassed) = obj.mu_0;
                %%% Reset sliding velocity
                obj.vs(hasPassed) = 0;
                %%% Reset sliding condition %%%
                obj.slide(hasPassed) = false;
                %%% Reset deformation velocities %%%
                obj.vx(hasPassed) = 0;
                obj.vy(hasPassed) = 0;
                %%% Reset relative velocities %%%
                obj.vrx(hasPassed) = 0;
                obj.vry(hasPassed) = 0;

                obj.theta_2(hasPassed) = -pi;           % Angle between sliding velocity and horizontal
            end
        end

        function obj = calculateVertical(obj, z_disp, z_dot)
            arguments
                obj 
                z_disp (:, :) single
                z_dot  (:, :) single
            end

            Press_ind = obj.hasPress;
            Z_eval = max(0, z_disp(Press_ind));     % Road cannot displace tyre downward
            Z_dot_eval = max(0, z_dot(Press_ind));  % Road cannot pull tyre downward

            pressVal = obj.kz(Press_ind) .* Z_eval + obj.cz .* Z_dot_eval;
            obj.press(Press_ind) = obj.initPress(Press_ind) + pressVal;
        end

        function obj = integrateDynamics(obj, dt, t, X_vec, method_name, omega, omega_z, re, v0, alpha)
        %INTEGRATEDYNAMICS is a simple wrapper function that takes in the user's
        %pre-defined dynamics as a function handle and integrates it
        %forward in time using the method specified in the method_name input
            % INPUTS
            % func:         Function handle of user's dynamics. Assumed to be
            %               only a function of (t, X).
            % dt:           Time step used. Will return next value at t + dt
            % t:            Current time value
            % X_vec:        Current state vector
            % method_name:  String containing the method of choice
            %
            % OUTPUTS
            % X_next:       The integrated states at time (t + dt)
            %
        
            arguments
                obj            
                dt              (1,1) single
                t               (1,1) single
                X_vec           (:,1) single
                method_name     
                omega   
                omega_z 
                re      
                v0      
                alpha   
            end
        
            switch lower(method_name)
                case 'euler'
                    obj.evaluateEuler_Brush(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                case 'implicit_euler'
                    obj = obj.evaluateImplicitEuler(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                % case 'adaptive_heun'
                %     t_current = t;
                %     t_target = t + dt;
                %     h_current = dt;
                %     while t_current < t_target
                %         % Calculate required step
                %         h_current = min(h_current,  t_target - t_current);
                % 
                %         % Call the adaptive step function
                %         [X_next, h_next] = adaptiveHeun_Brush(func, h_current, t_current, X_vec, brush_obj);
                % 
                %         % Update current time based on the step taken
                %         t_current = t_current + h_current;
                %         % Update time step
                %         h_current = h_next;
                %         % Update solution
                %         X_vec = X_next;
                %     end
                case 'verlet'
                    obj = obj.evaluateVerlet(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                case 'velocity_verlet'
                    obj = obj.evaluateVelocityVerlet(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                case 'tr_bdf2'
                    obj = obj.evaluateTRBDF2(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                case 'rk4'
                    obj.evaluateRK4_Brush(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                case 'rkf5'
                    obj = obj.evaluateRKF5_brush(dt, t, X_vec, omega, omega_z, re, v0, alpha);
                % case 'adaptive_rk45'
                %     t_current = t;
                %     t_target = t + dt;
                %     h_current = dt;
                %     for i = 1:1000
                %         if t_current < t_target
                %             % Calculate required step
                %             h_current = min(h_current,  t_target - t_current);
                % 
                %             % Call the adaptive step function
                %             [X_tp1, h_next] = adaptiveRK45(func, h_current, t_current, X_vec);
                % 
                %             % Update current time based on the step taken
                %             t_current = t_current + h_current;
                %             % Update time step
                %             h_current = h_next;
                %             % Update solution
                %             X_vec = X_tp1;
                %         end
                %     end
                otherwise
                    error('Unrecognised integration method: %s', method_name);
            end

            % Evaluate dynamical system at new state
            obj = obj.evaluateDynamics(t, omega, omega_z, re, v0, alpha);
        
        end

        function obj = evaluateEuler_Brush(obj, dt, t, X_vec, omega, omega_z, re, v0, alpha)
            k1= BrushVec_CPP.brushDynamics(t, X_vec, obj, omega, omega_z, re, v0, alpha);

            obj.X_next = X_vec + dt * k1;
            
        end

        function obj = evaluateRK4_Brush(obj, dt, t, X, omega, omega_z, re, v0, alpha)

            % [k1] = func(t, X, slide_obj, args);
            % [k2] = func(t + 0.5 *dt, X + 0.5 * dt * k1, slide_obj, args);
            % [k3] = func(t + 0.5 * dt, X + 0.5 * dt * k2, slide_obj, args);
            % [k4] = func(t + dt, X + dt * k3, slide_obj, args);
            % % Calculate the next state vector
            % X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
            
            % --- RK5(4)7S Method ---
            % Coefficients for the k_i stages
            A = [0; 2/9; 1/3; 5/9; 2/3; 1; 1];
            C = [0     , 0      , 0    , 0       , 0    , 0;
                 2/9   , 0      , 0    , 0       , 0    , 0;
                 1/12  , 1/4    , 0    , 0       , 0    , 0;
                 55/324, -25/108, 50/81, 0       , 0    , 0;
                 83/330, -13/22 , 61/66, 9/110   , 0    , 0;
                 -19/28, 9/4    , 1/7  , -27/7   , 22/7 , 0;
                 19/200, 0      , 3/5  , -243/400, 33/40, 7/80];
        
            % Coefficients for the 4th and 5th order solutions
            B_4 = [19/200  , 0, 3/5    , -243/400   , 33/40   , 7/80    , 0    ];
            % B_5 = [431/5000, 0, 333/500, -7857/10000, 957/1000, 193/2000, -1/50 ];
            % Calculate all k_i values
            k1 = BrushVec_CPP.brushDynamics(t + A(1)*dt, X, obj, omega, omega_z, re, v0, alpha);
            k2 = BrushVec_CPP.brushDynamics(t + A(2)*dt, X + C(2,1)* dt * k1, obj, omega, omega_z, re, v0, alpha);
            k3 = BrushVec_CPP.brushDynamics(t + A(3)*dt, X + C(3,1)* dt * k1 + C(3,2)*dt * k2, obj, omega, omega_z, re, v0, alpha);
            k4 = BrushVec_CPP.brushDynamics(t + A(4)*dt, X + C(4,1)* dt * k1 + C(4,2)*dt * k2 + C(4,3)* dt * k3, obj, omega, omega_z, re, v0, alpha);
            k5 = BrushVec_CPP.brushDynamics(t + A(5)*dt, X + C(5,1)* dt * k1 + C(5,2)*dt * k2 + C(5,3)* dt * k3 + C(5,4)* dt * k4, obj, omega, omega_z, re, v0, alpha);
            k6 = BrushVec_CPP.brushDynamics(t + A(6)*dt, X + C(6,1)* dt * k1 + C(6,2)*dt * k2 + C(6,3)* dt * k3 + C(6,4)* dt * k4 + C(6,5)* dt * k5, obj, omega, omega_z, re, v0, alpha);
            % k7 = func(t + A(7)*dt, X + C(7,1)* dt * k1 + C(7,2)*dt * k2 + C(7,3)* dt * k3 + C(7,4)* dt * k4 + C(7,5)* dt * k5 + C(7,6)* dt * k6, slide_obj, args);
            % Calculate the next state vector
            obj.X_next = X + B_4(1)* dt * k1 + B_4(2)* dt * k2 + B_4(3)* dt * k3 + B_4(4)* dt * k4 + B_4(5)* dt * k5 + B_4(6)* dt * k6;% + B_4(7)* dt * k7;
        
            % % % Update slide_obj with the final, correct state
            % % % Should call the calculateStresses method one last time
            % % % % [~, updated_obj] = func(t + dt, X_next, slide_obj);
            % % updated_obj = obj_k4;
        end

        function obj = updateStates(obj)
            coder.noImplicitExpansionInFunction();
            % Extract states into variables
            maxDelta_X = obj.mu(1:obj.numBrushes, 1) .* obj.press(1:obj.numBrushes) ./ obj.kx;
            maxDelta_Y = obj.mu(1:obj.numBrushes, 1) .* obj.press(1:obj.numBrushes) ./ obj.ky;
            
            obj.delta_x(1:obj.numBrushes, 1) = min( obj.X_next(1:obj.numBrushes, 1), maxDelta_X);
            obj.delta_y(1:obj.numBrushes, 1) = min( obj.X_next(obj.numBrushes + 1:2 * obj.numBrushes, 1), maxDelta_Y);
            obj.vx(1:obj.numBrushes, 1)      = obj.X_next(2 * obj.numBrushes + 1:3 * obj.numBrushes, 1);
            obj.vy(1:obj.numBrushes, 1)      = obj.X_next(3 * obj.numBrushes + 1:4 * obj.numBrushes, 1);
        end

        function obj = evaluateDynamics(obj, t, omega, omega_z, re, v0, alpha)
            
            obj = obj.updateStates();
       
            % Extract indices where pressure is being applied
            hasPressInd = obj.hasPress;
            slidInd = obj.slide;          
            % Only consider indices which has pressure applied
            isSliding = logical(slidInd .* hasPressInd);
            isSticking = logical((~slidInd) .* hasPressInd);

            % - calculate stresses at current time step based on difference between deformation velocity vx and sliding velocity vs_x
            % ---------------------------------------------------------------------
            %               Calculation for all elements with pressure
            % ---------------------------------------------------------------------
            % 	- calculate vrx and vry from algebraic constraints
            obj.vrx(hasPressInd) = omega(t) .* re + omega_z .* (obj.y(hasPressInd) + obj.delta_y(hasPressInd)) - 1e3 * v0(t) .* cos(alpha);
            obj.vry(hasPressInd) = -omega_z .* (obj.x(hasPressInd) + obj.delta_x(hasPressInd)) - 1e3 * v0(t) .* sin(alpha);
       
            % ---------------------------------------------------------------------
            %               Sliding Region
            % ---------------------------------------------------------------------
            % 	if sliding and has pressure
            % 	- calculate vs_x, vs_y in order to calculate vs (vs_x = vx - vrx/ vs_y = vy - vry) [ ONLY FOR ELEMENTS THAT ARE SLIDING ]
            obj.vs_x(isSliding) = obj.vx(isSliding) - obj.vrx(isSliding);
            obj.vs_y(isSliding) = obj.vy(isSliding) - obj.vry(isSliding);
            obj.vs(isSliding) = hypot(obj.vs_x(isSliding), obj.vs_y(isSliding));
            % 	- calculate theta_2 = atan2(vs_y, vs_x)
            obj.theta_2(isSliding) = pi + atan2(obj.vs_y(isSliding), obj.vs_x(isSliding)); % Minus pi to constrain it to collinear with velocity
            % 	else
            % ---------------------------------------------------------------------
            %               Adhesion Region
            % ---------------------------------------------------------------------
            % 	- vs_x, vs_y = 0
            obj.vs_x(isSticking) = 0;
            obj.vs_y(isSticking) = 0;
            obj.vs(isSticking) = 0;
            % 	- vx, vy = vrx, vry
            obj.vx(isSticking) = obj.vrx(isSticking);
            obj.vy(isSticking) = obj.vry(isSticking);
            % 	- calculate theta_2 = atan2(vry, vrx)
            obj.theta_2(isSticking) = atan2(obj.vry(isSticking), obj.vrx(isSticking)); % Minus pi to constrain it to collinear with velocity
            % 	endif
            % ---------------------------------------------------------------------
            %               Calculation for all elements with pressure
            % ---------------------------------------------------------------------
            % 	- calculate friction coefficient based on vs and press from mastercurve equation
            obj.mu(hasPressInd) = BrushVec_CPP.compute_friction_coefficient(obj.press(hasPressInd), obj.p_0, obj.p_ref, obj.q, ...
                                                                         obj.mu_0,obj.mu_m, obj.h,...
                                                                         obj.vs(hasPressInd), obj.v_m);
            % 	- calculate stress: tauX = mu * press * cos(theta_2) / tauY = mu * press * sin(theta_2);
            obj.tauX(hasPressInd) = obj.mu(hasPressInd) .* obj.press(hasPressInd) .* cos(obj.theta_2(hasPressInd));
            obj.tauY(hasPressInd) = obj.mu(hasPressInd) .* obj.press(hasPressInd) .* sin(obj.theta_2(hasPressInd));

            % Set the values that has no pressure to zero
            obj.vx(~hasPressInd) = 0;
            obj.vy(~hasPressInd) = 0;
            obj.vs_x(~hasPressInd) = 0;
            obj.vs_y(~hasPressInd) = 0;
            obj.vs(~hasPressInd) = 0;
            obj.tauX(~hasPressInd) = 0;
            obj.tauY(~hasPressInd) = 0;
      
        end
    end


    methods (Static)

        function [dX] = brushDynamics(t, X, obj, omega, omega_z, re, v0, alpha)
        %BRUSHDYNAMICS Summary of this function goes here
            % Get number of brushes in vector
            num_masses = obj.numBrushes;
            % Extract states into variables
            DeltaX = X(1:num_masses);
            DeltaY = X(num_masses + 1:2 * num_masses);
            VX     = X(2 * num_masses + 1:3 * num_masses);
            VY     = X(3 * num_masses + 1:4 * num_masses);
       
            % Extract indices where pressure is being applied
            hasPress = obj.hasPress;
            slidInd = obj.slide;
            % Load variables into temperary variables
            
            vrx     = obj.vrx    ; 
            vry     = obj.vry    ; 
            vs_x    = obj.vs_x   ; 
            vs_y    = obj.vs_y   ; 
            vs      = obj.vs     ;
            tauX    = obj.tauX   ; 
            tauY    = obj.tauY   ; 
            mu      = obj.mu     ;
            theta_2 = obj.theta_2;
            % Only consider indices which has pressure applied
            isSliding = logical(slidInd .* hasPress);
            isSticking = logical((~slidInd) .* hasPress);

            % - calculate stresses at current time step based on difference between deformation velocity vx and sliding velocity vs_x
            % ---------------------------------------------------------------------
            %               Calculation for all elements with pressure
            % ---------------------------------------------------------------------
            % 	- calculate vrx and vry from algebraic constraints
            vrx(hasPress) = omega(t) .* re + omega_z .* (obj.y(hasPress) + DeltaX(hasPress)) - 1e3 * v0(t) .* cos(alpha);
            vry(hasPress) = -omega_z .* (obj.x(hasPress) + DeltaY(hasPress)) - 1e3 * v0(t) .* sin(alpha);
       
            % ---------------------------------------------------------------------
            %               Sliding Region
            % ---------------------------------------------------------------------
            % 	if sliding and has pressure
            % 	- calculate vs_x, vs_y in order to calculate vs (vs_x = vx - vrx/ vs_y = vy - vry) [ ONLY FOR ELEMENTS THAT ARE SLIDING ]
            vs_x(isSliding) = VX(isSliding) - vrx(isSliding);
            vs_y(isSliding) = VY(isSliding) - vry(isSliding);
            vs(isSliding) = hypot(vs_x(isSliding), vs_y(isSliding));
            % 	- calculate theta_2 = atan2(vs_y, vs_x)
            theta_2(isSliding) = atan2(vs_y(isSliding), vs_x(isSliding)) - pi; % Minus pi to constrain it to collinear with velocity
            % 	else
            % ---------------------------------------------------------------------
            %               Adhesion Region
            % ---------------------------------------------------------------------
            % 	- vs_x, vs_y = 0
            vs(isSticking) = 0;
            % 	- vx, vy = vrx, vry
            VX(isSticking) = vrx(isSticking);
            VY(isSticking) = vry(isSticking);
            % 	- calculate theta_2 = atan2(vry, vrx)
            theta_2(isSticking) = atan2(vry(isSticking), vrx(isSticking)) + pi; % Minus pi to constrain it to collinear with velocity
            % 	endif
            % ---------------------------------------------------------------------
            %               Calculation for all elements with pressure
            % ---------------------------------------------------------------------
            % 	- calculate friction coefficient based on vs and press from mastercurve equation
            mu(hasPress) = BrushVec_CPP.compute_friction_coefficient(obj.press(hasPress), obj.p_0, obj.p_ref, obj.q, ...
                                                                         obj.mu_0,obj.mu_m, obj.h,...
                                                                         vs(hasPress), obj.v_m);
            mu = mu(:);
            % 	- calculate stress: tauX = mu * press * cos(theta_2) / tauY = mu * press * sin(theta_2);
            tauX(hasPress) = mu(hasPress) .* obj.press(hasPress) .* cos(theta_2(hasPress));
            tauY(hasPress) = mu(hasPress) .* obj.press(hasPress) .* sin(theta_2(hasPress));

            % Set the values that has no pressure to zero
            VX(~hasPress) = 0;
            VY(~hasPress) = 0;
            tauX(~hasPress) = 0;
            tauY(~hasPress) = 0;
            % - construct derivative of state vector: dX = [vx; (tauX - k * delta_x - c * vx)/m; vy; (tauY - k * delta_y - c * vy)/m];
            dVx = (tauX - obj.kx * DeltaX - obj.cx * VX);
            dVy = (tauY - obj.ky * DeltaY - obj.cy * VY);
            dVx(isSticking) = 0;
            dVy(isSticking) = 0;
            
            dX = [VX;
                  dVx;
                  VY;
                  dVy;
                  ];
              
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
                tauX   (:,1)  single
                tauY   (:,1)  single
                mu_0   (1,1)    single
                press  (:,1)  single
            end
            
            tau = hypot(tauX, tauY);
            slide = (tau > mu_0 .* press);
        end

        function [hasPress] = has_pressure(press, p_0)
            % Determines if elements have vertical pressure applied
            %
            % Parameters:
            %   press:  Vertical Pressure
            %   p_0:    Minimum pressure value   
            %
            % Returns:
            %   hasPress - Logical array indicating elements with vertical
            %              pressure
            
            arguments
                press   (:,1)  single
                p_0     (1,1)  single
            end
            
            hasPress = (press > p_0);
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
                currentX    (:,1) single
                minX        (:,1) single
            end
            passed = currentX <= minX;
        end

    end
end

            