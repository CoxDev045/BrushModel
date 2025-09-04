function [dX, sticking_obj] = stickingDynamics(t,X, sticking_obj, forcing_args)
    % INPUTS
    % X:        Current state vector containing the following:
    %           [delta_x(:), delta_y(:), vx(:), vy(:)];
    % obj:      Brush object containing all the necessary
    %           properties for integration
    % args:     
    % OUTPUTS
    % dx:       Vector containing the derivatives of the next state
    %           [ddelta_x;
    %            ddelta_y;
    %            dvx;
    %            dvy];
    %    

    % Extract states from state vector
    num_masses = length(sticking_obj.vrx(~sticking_obj.slide));
    DeltaX = X(1:num_masses);
    DeltaY = X(num_masses + 1:2 * num_masses);
    VrX     = X(2 * num_masses + 1:3 * num_masses);
    VrY     = X(3 * num_masses + 1:4 * num_masses);

   % Update velocity history
    sticking_obj.prev3_vrx(~sticking_obj.slide)   = sticking_obj.prev2_vrx(~sticking_obj.slide);
    sticking_obj.prev3_vry(~sticking_obj.slide)   = sticking_obj.prev2_vry(~sticking_obj.slide);
    sticking_obj.prev2_vrx(~sticking_obj.slide)   = sticking_obj.prev1_vrx(~sticking_obj.slide);
    sticking_obj.prev2_vry(~sticking_obj.slide)   = sticking_obj.prev1_vry(~sticking_obj.slide);
    sticking_obj.prev1_vrx(~sticking_obj.slide)   = VrX;
    sticking_obj.prev1_vry(~sticking_obj.slide)   = VrY;

    % Extract additional parameters
    omega   = forcing_args.omega;
    re      = forcing_args.re;
    omega_z = forcing_args.omega_z;
    v0      = forcing_args.v0;
    alpha   = forcing_args.alpha;
    dt      = forcing_args.dt;

    % Calculate New Velocities
    sticking_obj.vrx(~sticking_obj.slide) = omega .* re + omega_z .* (sticking_obj.y(~sticking_obj.slide) + DeltaY) - v0 .* cos(alpha);
    sticking_obj.vry(~sticking_obj.slide) = -omega_z .* (sticking_obj.x(~sticking_obj.slide) + DeltaX) - v0 .* sin(alpha);

    % 2nd order backward difference for acceleration estimation
    sticking_obj.dvrx(~sticking_obj.slide) = (3 .* sticking_obj.vrx(~sticking_obj.slide) - 4 .* sticking_obj.prev1_vrx(~sticking_obj.slide) + sticking_obj.prev2_vrx(~sticking_obj.slide)) ./ (2 * dt);
    sticking_obj.dvry(~sticking_obj.slide) = (3 .* sticking_obj.vry(~sticking_obj.slide) - 4 .* sticking_obj.prev1_vry(~sticking_obj.slide) + sticking_obj.prev2_vry(~sticking_obj.slide)) ./ (2 * dt);
    
    
    % Construct derivative of state vector
    dX = [sticking_obj.vrx(~sticking_obj.slide);
          sticking_obj.vry(~sticking_obj.slide);
          sticking_obj.dvrx(~sticking_obj.slide);
          sticking_obj.dvry(~sticking_obj.slide);]; 
end
