function [dX,obj] = brushDynamics(t, X, obj, omega, omega_z, re, v0, alpha)
%BRUSHDYNAMICS Summary of this function goes here
    num_masses = obj.numBrushes;
    DeltaX = X(1:num_masses);
    DeltaY = X(num_masses + 1:2 * num_masses);
    VX     = X(2 * num_masses + 1:3 * num_masses);
    VY     = X(3 * num_masses + 1:4 * num_masses);
    % - calculate stresses at current time step based on difference between deformation velocity vx and sliding velocity vs_x
    % ---------------------------------------------------------------------
    %               Calculation for all elements
    % ---------------------------------------------------------------------
    % 	- calculate vrx and vry from algebraic constraints
    obj.vrx = omega .* re + omega_z .* (obj.y + DeltaX) - v0 .* cos(alpha);
    obj.vry = -omega_z .* (obj.x + DeltaY) - v0 .* sin(alpha);

    slidInd = obj.slide;
    % ---------------------------------------------------------------------
    %               Sliding Region
    % ---------------------------------------------------------------------
    % 	if sliding
    % 	- calculate vs_x, vs_y in order to calculate vs (vs_x = vx - vrx/ vs_y = vy - vry) [ ONLY FOR ELEMENTS THAT ARE SLIDING ]
    obj.vs_x(slidInd) = VX(slidInd) - obj.vrx(slidInd);
    obj.vs_y(slidInd) = VY(slidInd) - obj.vry(slidInd);
    obj.vs(slidInd) = hypot(obj.vs_x(slidInd), obj.vs_y(slidInd));
    % 	- calculate theta_2 = atan2(vs_y, vs_x)
    obj.theta_2(slidInd) = atan2(obj.vs_y(slidInd), obj.vs_x(slidInd)) - pi; % Minus pi to constrain it to collinear with velocity
    % 	else
    % ---------------------------------------------------------------------
    %               Adhesion Region
    % ---------------------------------------------------------------------
    % 	- vs_x, vs_y = 0
    obj.vs_x(~slidInd) = 0;
    obj.vs_y(~slidInd) = 0;
    obj.vs(~slidInd) = 0;
    % 	- vx, vy = vrx, vry
    VX(~slidInd) = obj.vrx(~slidInd);
    VY(~slidInd) = obj.vry(~slidInd);
    % 	- calculate theta_2 = atan2(vry, vrx)
    obj.theta_2(~slidInd) = atan2(obj.vry(~slidInd), obj.vrx(~slidInd)) - pi; % Minus pi to constrain it to collinear with velocity
    % 	endif
    % ---------------------------------------------------------------------
    %               Calculation for all elements
    % ---------------------------------------------------------------------
    % 	- calculate friction coefficient based on vs and press from mastercurve equation
    obj.mu = BrushVec_CPP.compute_friction_coefficient(obj.press, obj.p_0, obj.p_ref, obj.q, ...
                                                        obj.mu_0,obj.mu_m, obj.h, obj.vs, obj.v_m);
    % 	- calculate stress: tauX = mu * press * cos(theta_2) / tauY = mu * press * sin(theta_2);
    obj.tauX = obj.mu .* obj.press .* cos(obj.theta_2);
    obj.tauY = obj.mu .* obj.press .* sin(obj.theta_2);
    
    % - construct derivative of state vector: dX = [vx; (tauX - k * delta_x - c * vx)/m; vy; (tauY - k * delta_y - c * vy)/m];
    dX = [VX;
          (obj.tauX - obj.kx * obj.delta_x - obj.cx * VX);
          VY;
          (obj.tauY - obj.ky * obj.delta_y - obj.cy * VY);
          ];
      
end

