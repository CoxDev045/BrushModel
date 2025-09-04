function [dX, slide_obj] = slidingDynamics(t,X, slide_obj, args)
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
    num_masses = length(slide_obj.vrx(slide_obj.slide));
    DeltaX = X(1:num_masses);
    DeltaY = X(num_masses + 1:2 * num_masses);
    VX     = X(2 * num_masses + 1:3 * num_masses);
    VY     = X(3 * num_masses + 1:4 * num_masses);

    % Get forcing at time t
    slide_obj = slide_obj.calculateStresses(t, X, args);

    Fx = slide_obj.tauX(slide_obj.slide);
    Fy = slide_obj.tauY(slide_obj.slide);
    
    % Apply dynamics only where sliding occurs
    d_delta_x = VX;
    d_vx      = (Fx - slide_obj.kx .* DeltaX - slide_obj.cx .* VX) ./ slide_obj.m;

    d_delta_y = VY;
    d_vy      = (Fy - slide_obj.ky .* DeltaY - slide_obj.cy .* VY) ./ slide_obj.m;

    % Reassemble the output into a single column vector
    dX = [d_delta_x; d_delta_y; d_vx; d_vy];  
end
