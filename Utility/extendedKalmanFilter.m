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
    [numStates, ~] = size(Q);
    [numMeasured, ~] = size(R);
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

    if ~isempty(Y)
        % ---------- Actualisation Step --------------
        % Measurement innovation
        measuredState = measure(t, x_pred, numStates, numMeasured);
        V_pred = Y - measuredState;
    
        % Covariance innovation
        % Calculate jacobian based on predicted state
        [H_pred,facH] = numjac(@(T, X) measure(T, X, numStates, numMeasured) , t,x_pred,measuredState,thresh_scal, facH);
        S_pred = H_pred * P_pred * H_pred.' + R; % Constant measurement noise covariance
    
        % Kalman Gain
        % KG * S = P * H.'
        % KG = (P * H.') * inv(S)

        % S_pred: size N x N
        b = P_pred * H_pred.'; % size N x M (N states by M measurements)
        
        % C = cond(S_pred);
        if cond(S_pred) < 1e12
            % Slightly faster and can be more accurate than inv(S)
            S_decomp = decomposition(S_pred);
            % Perform gauss elimination 
            KG = b / S_decomp;
        else
            % Only use psuedo inverse when condition number much greater
            % than 1
            % Psuedo inverse is costly to compute and should not be default
            KG = b * pinv(S_pred);
            % --- Alternative ----
            % Use least squares to find inverse of S_pred that solves the
            % system x * A = b
            % KG = lsqminnorm(S_pred, P_pred * H_pred.');
        
        end
        % Actualisation of the state estimate
        X_next = x_pred + KG * V_pred;
        % % ------- Constrain parameters to search space -----
        % % Constrain the parameters to [1e-8, 1]
        % X_next(3:end) = max(1e-8, min(20, X_next(3:end)));
    
        % Actualisation of the covariance estimate
        PC = (eye(size(P_pred)) - KG * H_pred) * P_pred;
    else
        X_next = x_pred;
        % % ------- Constrain parameters to search space -----
        % % Constrain the parameters to [1e-8, 1]
        % X_next(3:end) = max(1e-8, min(1, X_next(3:end)));
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
    
end


function measuredState = measure(t, X, numStates, numMeasured)
    % Measure the first state for now...
    % h = [1, 0, 0, 0, 0;
    %      0, 1, 0, 0, 0;];
    h = eye(numMeasured, numStates);
    % h = fliplr(h); % Measure only vel
    measuredState = h * X;
end