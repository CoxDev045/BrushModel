function [X_next] = evaluateRK4_Brush(func, dt, t, X, slide_obj, args)

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
    k1 = func(t + A(1)*dt, X, slide_obj, args);
    k2 = func(t + A(2)*dt, X + C(2,1)* dt * k1, slide_obj, args);
    k3 = func(t + A(3)*dt, X + C(3,1)* dt * k1 + C(3,2)*dt * k2, slide_obj, args);
    k4 = func(t + A(4)*dt, X + C(4,1)* dt * k1 + C(4,2)*dt * k2 + C(4,3)* dt * k3, slide_obj, args);
    k5 = func(t + A(5)*dt, X + C(5,1)* dt * k1 + C(5,2)*dt * k2 + C(5,3)* dt * k3 + C(5,4)* dt * k4, slide_obj, args);
    k6 = func(t + A(6)*dt, X + C(6,1)* dt * k1 + C(6,2)*dt * k2 + C(6,3)* dt * k3 + C(6,4)* dt * k4 + C(6,5)* dt * k5, slide_obj, args);
    % k7 = func(t + A(7)*dt, X + C(7,1)* dt * k1 + C(7,2)*dt * k2 + C(7,3)* dt * k3 + C(7,4)* dt * k4 + C(7,5)* dt * k5 + C(7,6)* dt * k6, slide_obj, args);
    % Calculate the next state vector
    X_next = X + B_4(1)* dt * k1 + B_4(2)* dt * k2 + B_4(3)* dt * k3 + B_4(4)* dt * k4 + B_4(5)* dt * k5 + B_4(6)* dt * k6;% + B_4(7)* dt * k7;

    % % % Update slide_obj with the final, correct state
    % % % Should call the calculateStresses method one last time
    % % % % [~, updated_obj] = func(t + dt, X_next, slide_obj);
    % % updated_obj = obj_k4;
end