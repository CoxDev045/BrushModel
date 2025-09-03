function [y_next, h_next] = adaptiveRK45(func, dt, t_current, X_vec)
    % adaptive_ODE performs one step of numerical integration.
    % An adaptive timestep with RK4 as first step and RK5 as
    % finer resolution step.
    %
    %   f: Function handle for the ODE: dy/dt = f(t, y)
    %   dt: Current time step
    %   t_current: Current time
    %   t_target: Final time step to solve for
    %   X_vec: Current solution vector
    %   args: Arguments for function func(t, x, args)
    %
    %   y_next: Solution at time t+h
    %   h_next: Time step at the end of the current step 

    % % % --- RK45-Sarafyan Method ---
    % % % Coefficients for the k_i stages
    % % A = [0; 1/2; 1/2; 1; 2/3; 2/10];
    % % C = [0     , 0       , 0      , 0     , 0;
    % %      1/2   , 0       , 0      , 0     , 0;
    % %      1/4   , 1/4     , 0      , 0     , 0;
    % %      0     , -1      , 2      , 0     , 0;
    % %      7/27  , 10/27   , 0      , 1/27  , 0;
    % %      28/625, -125/625, 546/625, 54/625, -378/625];
    % % 
    % % % Coefficients for the 4th and 5th order solutions
    % % B_4 = [1/6, 0, 4/6, 3, 0, 0];
    % % B_5 = [14/336, 0, 0, 35/336, 162/336, 125/336];

    % % % --- RK5(4)7M Method ---
    % % % Coefficients for the k_i stages
    % % A = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
    % % C = [0         , 0          , 0         , 0       , 0          , 0;
    % %      1/5       , 0          , 0         , 0       , 0          , 0;
    % %      3/40      , 9/40       , 0         , 0       , 0          , 0;
    % %      44/45     , -56/15     , 32/9      , 0       , 0          , 0;
    % %      19372/6561, -25360/2187, 64448/6561, -212/729, 0          , 0;
    % %      -9017/3168, -355/33    , 46732/5247, 49/176  , -5103/18656, 0;
    % %      35/384    , 0          , 500/1113  , 125/192 , -2187/6784 , 11/84];
    % % 
    % % % Coefficients for the 4th and 5th order solutions
    % % B_4 = [35/384    , 0, 500/1113  , 125/192, -2187/6784   , 11/84   , 0    ];
    % % B_5 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40 ];

    % % % --- RK5(4)7S Method ---
    % % % Coefficients for the k_i stages
    % % A = [0; 2/9; 1/3; 5/9; 2/3; 1; 1];
    % % C = [0     , 0      , 0    , 0       , 0    , 0;
    % %      2/9   , 0      , 0    , 0       , 0    , 0;
    % %      1/12  , 1/4    , 0    , 0       , 0    , 0;
    % %      55/324, -25/108, 50/81, 0       , 0    , 0;
    % %      83/330, -13/22 , 61/66, 9/110   , 0    , 0;
    % %      -19/28, 9/4    , 1/7  , -27/7   , 22/7 , 0;
    % %      19/200, 0      , 3/5  , -243/400, 33/40, 7/80];
    % % 
    % % % Coefficients for the 4th and 5th order solutions
    % % B_4 = [19/200  , 0, 3/5    , -243/400   , 33/40   , 7/80    , 0    ];
    % % B_5 = [431/5000, 0, 333/500, -7857/10000, 957/1000, 193/2000, -1/50 ];

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


    % Initialise h and y
    y = X_vec;
    h_current = dt;
    
    % Adaptive time step parameters
    rtol = 1e-14;
    atol = 1e-13;
    h_min = sqrt(eps);
    h_max = 1e-3;
    max_retries = 20;
    tolerance = max(rtol * norm(y), atol);

    retries = 0;
    while retries < max_retries
        % Calculate all k_i values
        k1 = h_current * func(t_current + A(1)*h_current, y);
        k2 = h_current * func(t_current + A(2)*h_current, y + C(2,1)*k1);
        k3 = h_current * func(t_current + A(3)*h_current, y + C(3,1)*k1 + C(3,2)*k2);
        k4 = h_current * func(t_current + A(4)*h_current, y + C(4,1)*k1 + C(4,2)*k2 + C(4,3)*k3);
        k5 = h_current * func(t_current + A(5)*h_current, y + C(5,1)*k1 + C(5,2)*k2 + C(5,3)*k3 + C(5,4)*k4);
        k6 = h_current * func(t_current + A(6)*h_current, y + C(6,1)*k1 + C(6,2)*k2 + C(6,3)*k3 + C(6,4)*k4 + C(6,5)*k5);
        % k7 = h_current * func(t_current + A(7)*h_current, y + C(7,1)*k1 + C(7,2)*k2 + C(7,3)*k3 + C(7,4)*k4 + C(7,5)*k5 + C(7,6)*k6, args{:});
        % Calculate the two solutions (4th and 5th order)
        y_RK4 = y + B_4(1)*k1 + B_4(2)*k2 + B_4(3)*k3 + B_4(4)*k4 + B_4(5)*k5 + B_4(6)*k6;
        y_RK5 = y + B_5(1)*k1 + B_5(2)*k2 + B_5(3)*k3 + B_5(4)*k4 + B_5(5)*k5 + B_5(6)*k6;

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
    y_next = y; % Return the original point
    h_next = h_current;
end
