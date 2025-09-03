function [X_next, KG, PC, output] = extendedKalmanFilter(func, t, X, Y, options)
%EXTENDEDKALMANFILTER uses the EKF formulation to step a non-linear
%function forward by 1 time step
%   func: Model dynamics function. Must output a column vector of states.
%   t: current time
%   X: Current state vector
%   Y: Current Measurment vector
%   options:    ThreshScal: Threshold scaling for numjac
%               facF: Used for adaptive step sizing in numjac (For dynamics)
%               facH: Used for adaptive step sizing in numjac (For measurements)
%               ProcessNoiseCov: Process noise covariance
%               MeasurementNoiseCov: Measurement noise covariance
%               ErrorCov: Error covariance. Used to quantify error
%               dt: Current timestep
%               ControlVec: Control vector (used for propogation without measurements)
%
%   X_next: Predicted state vector
%   KG: Kalman gain matrix
%   PC: Prediction Covariance matrix
    
    % Extract options
    thresh_scal = options.ThreshScal;
    facF = options.facF;
    facH = options.facH;
    Q = options.ProcessNoiseCov;
    R = options.MeasurementNoiseCov;
    P = options.ErrorCov;
    dt = options.dt;
    u = options.ControlVec;
    % ------------- Prediction Step ---------------
    % Predict state based on previous state and dynamics model
    % Integrate dynamical system to give predicted state
    % func is assumed to handle dynamics and integration and output system
    % states (e.g displacement and velocity)
    % Calculate jacobian based on predicted state
    [F_prev,facF] = numjac(func,t,X,X,thresh_scal, facF);
    F_xx_prev = F_prev.' * F_prev;

    x_pred = func(t, X) + 0.5 * trace(F_xx_prev * P);

    % Calculate jacobian based on predicted state
    [F_pred,facF] = numjac(func,t,X,x_pred,thresh_scal, facF);

    % Predict covariance estimate
    P_pred = F_pred * P * F_pred.' + ...
             0.5 * trace(F_xx_prev * P * F_xx_prev) + ...
             Q; % Constant process noise covariance

    if isempty(u)
        % ---------- Actualisation Step --------------
        % Measurement innovation
        measuredState = measure(t, x_pred);
        % Calculate jacobian based on predicted state
        [H_pred,facH] = numjac(@measure,t,x_pred,measuredState,thresh_scal, facH);
        H_xx_pred = H_pred.' * H_pred;
        % Error term
        V_pred = Y - measuredState - 0.5 * trace(H_xx_pred * P_pred);
    
        % Covariance innovation
        S_pred = H_pred * P_pred * H_pred.' + ...
                 0.5 * trace(H_xx_pred * P_pred * H_xx_pred.') + ...
                 R; % Constant measurement noise covariance
    
        % Kalman Gain
        b = P_pred * H_pred.'; % size N x M (N states by M measurements)
        % S_pred: size N x N
        % KG * S = P * H.'
        % ------- Solve the following system ------
        % KG = (P * H.') * inv(S)

        C = cond(S_pred);
        if C < 1e12
            dA = decomposition(S_pred);
            KG = b / dA;
        else
            KG = b * pinv(S_pred);
        end
        % KG = lsqminnorm(S_pred, P_pred * H_pred.');
        % Actualisation of the state estimate
        X_next = x_pred + KG * V_pred;
        % Constrain the parameters to [1e-8, 1]
        X_next(3:end) = max(1e-8, min(20, X_next(3:end)));
    
        % Actualisation of the covariance estimate
        PC = (eye(size(P_pred)) - KG * H_pred) * P_pred;
    else
        X_next = x_pred;
        % Constrain the parameters to [1e-8, 1]
        X_next(3:end) = max(1e-8, min(1, X_next(3:end)));
        PC = P_pred;
        KG = [];
    end
    % ---------- Save output for next round -----------
    output.facF = facF; % For adaptive step size in numjac (for F_pred)
    output.facH = facH; % For adaptive step size in numjac (for H_pred)
    output.ThreshScal = thresh_scal;
    output.ProcessNoiseCov = Q;
    output.MeasurementNoiseCov = R;
    output.ErrorCov = PC;
    output.dt = dt;
    output.ControlVec = u;
    
end


function measuredState = measure(t, X)
    % Measure the first state for now...
    h = [1, 0, 0, 0, 0;
         0, 1, 0, 0, 0;];
    measuredState = h * X;
end