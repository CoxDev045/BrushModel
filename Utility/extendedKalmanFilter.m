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
    x_pred = func(t, X);

    % Calculate jacobian based on predicted state
    [F_pred,facF] = numjac(func,t,X,x_pred,thresh_scal, facF);

    % Predict covariance estimate
    P_pred = F_pred * P * F_pred.' + Q; % Constant process noise covariance

    if isempty(u)
        % ---------- Actualisation Step --------------
        % Measurement innovation
        measuredState = measure(t, x_pred);
        V_pred = Y - measuredState;
    
        % Covariance innovation
        % Calculate jacobian based on predicted state
        [H_pred,facH] = numjac(@measure,t,x_pred,measuredState,thresh_scal, facH);
        S_pred = H_pred * P_pred * H_pred.' + R; % Constant measurement noise covariance
    
        % Kalman Gain
        S_decomp = chol(S_pred);
        KG = (P_pred * H_pred.' / S_decomp) / S_decomp.';
    
        % Actualisation of the state estimate
        X_next = x_pred + KG .* V_pred;
    
        % Actualisation of the covariance estimate
        PC = (eye(size(P_pred)) - KG * H_pred) * P_pred;
    else
        X_next = x_pred;
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