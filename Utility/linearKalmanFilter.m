function [X_next, KG, PC, output] = linearKalmanFilter(t, X, Y, options)
%EXTENDEDKALMANFILTER uses the EKF formulation to step a non-linear
%function forward by 1 time step
%   func: Model dynamics function. Must output a column vector of states.
%   t: current time
%   X: Current state vector
%   Y: Current Measurment vector
%   options:    StateTransition: State transition matrix
%               MeasurementModel: Measurement Matrix
%               ProcessNoiseCov: Process noise covariance
%               MeasurementNoiseCov: Measurement noise covariance
%               ErrorCov: Error covariance. Used to quantify error
%               dt: Current timestep
%
%   X_next: Predicted state vector
%   KG: Kalman gain matrix
%   PC: Prediction Covariance matrix
    
    % Extract options
    F_pred = options.StateTransition;
    B = options.InputMatrix;
    u = options.ControlMatrix;
    H_pred = options.MeasurementModel;
    Q = options.ProcessNoiseCov;
    R = options.MeasurementNoiseCov;
    P = options.ErrorCov;
    % ------------- Prediction Step ---------------
    % Predict state based on previous state and dynamics model
    % Integrate dynamical system to give predicted state
    % func is assumed to handle dynamics and integration and output system
    % states (e.g displacement and velocity)
    if ~isempty(u)
        x_pred = F_pred * X + B * u(t);
    else
        x_pred = F_pred * X;
    end
    % Predict covariance estimate
    P_pred = F_pred * P * F_pred.' + Q; % Constant process noise covariance

    if ~isempty(Y)
        % ---------- Actualisation Step --------------
        % Measurement innovation
        measuredState = H_pred * x_pred;
        V_pred = Y - measuredState;
    
        % Covariance innovation
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
    
        % Actualisation of the covariance estimate
        % Using Joseph form of Covariance (Better numerical Stability)
        % % PC = (eye(size(P_pred)) - KG * H_pred) * P_pred;
        % Reference: https://kalman-filter.com/joseph-form/
        KG_H = KG * H_pred;
        A = eye(size(KG_H)) - KG_H;
        PC = A * P_pred * A.' + KG * R * KG.';
    else
        X_next = x_pred;
        PC = P_pred;
        KG = [];
    end
    % ---------- Save output for next round -----------
    output.StateTransition = F_pred;
    output.ControlMatrix = B;
    output.Forcing = u;
    output.MeasurementModel = H_pred;
    output.ProcessNoiseCov = Q;
    output.MeasurementNoiseCov = R;
    output.ErrorCov = PC;
    
end
