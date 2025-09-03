function [X_next, KG, PC, output] = unscentedKalmanFilter(func, t, X, Y, options)
%UNSCENTEDKALMANFILTER uses the UKF formulation to step a non-linear
%function forward by 1 time step
%   func: Model dynamics function. Must output a column vector of states.
%   t: current time
%   X: Current state vector
%   Y: Current Measurment vector
%   options:    ProcessNoiseCov: Process noise covariance
%               MeasurementNoiseCov: Measurement noise covariance
%               ErrorCov: Error covariance. Used to quantify error
%               dt: Current timestep
%               ControlVec: Control vector (used for propogation without measurements)
%
%   X_next: Predicted state vector
%   KG: Kalman gain matrix
%   PC: Prediction Covariance matrix
    
    % Extract options
    Q = options.ProcessNoiseCov;
    R = options.MeasurementNoiseCov;
    P = options.ErrorCov;
    dt = options.dt;

    [numStates, ~] = size(Q);
    [numMeasured, ~] = size(R);

    n = numStates;
    % Predefined variables
    alpha = 1;
    kappa = 3 - n;
    beta = 2;
    % ------------- Calculate weights ---------------
    % Calculate lambda
    lambda = alpha^2 * (n + kappa) - n;
    % Calculate weights
    W_mean = zeros(2 * n +1, 1);
    W_cov = zeros(2 * n +1, 1); 
    W_mean(1) = lambda/(n + lambda);
    W_cov(1) = W_mean(1) + (1 - alpha^2 + beta);
    W_mean(2:end) = 1 / (2 * (n + lambda));
    W_cov(2:end) = 1 / (2 * (n + lambda));
    % ------------- Form Sigma Points ---------------
    X_sigma = formSigmaPoints(X, P, n, lambda);
    % ------------- Propogate through function ---------------
    X_pred = zeros(numStates, 2 * n + 1);
    for i = 1:2 * n + 1
        X_pred(:, i) = func(t, X_sigma(:, i));
    end
    % ------------- Prediction Step ---------------
    % Calculate weighted average and sum along the rows
    m_pred = sum(W_mean.' .* X_pred, 2);
    P_pred = zeros(size(P));
    for i = 1:2 * n + 1
        P_weighted = (W_cov(i) * (X_pred(:, i) - m_pred) * (X_pred(:, i) - m_pred).');
        P_pred = P_pred + P_weighted;
    end
    P_pred = P_pred + Q;

    if ~isempty(Y)
        % ------------- Form sigma points for update step ---------------
        X_sigma_pred = formSigmaPoints(m_pred, P_pred, n, lambda);
        
        % ---------- Actualisation Step --------------
        % Measurement innovation
        Z = zeros(numMeasured, 2 * n + 1);
        for i = 1:2 * n + 1
            Z(:, i) = measure(t, X_sigma_pred(:, i), numStates, numMeasured);
        end
        
        % Calculate mean of measured states
        mu_pred = sum(W_mean.' .* Z, 2);
        
        
        S_pred = zeros(numMeasured);
        C_pred = zeros(numStates, numMeasured);
        for i = 1:2 * n + 1
            % Calculate measurement covariance
            %                        (M x 1   - M x 1  ) * (M X 1   - M x 1)^T = M X M
            S_weighted = (W_cov(i) * (Z(:, i) - mu_pred) * (Z(:, i) - mu_pred).');
            S_pred = S_pred + S_weighted;
            % Calculate cross-covariance of measurement and state
            %                        (N x 1        - N x 1 ) * (M X 1   - M x 1)^T = N X M
            C_weighted = (W_cov(i) * (X_pred(:, i) - m_pred) * (Z(:, i) - mu_pred).');
            C_pred = C_pred + C_weighted;
        end
        S_pred = S_pred + R;
        
        % Compute Kalman Gain
        if cond(S_pred) < 1e12
            % Slightly faster and can be more accurate than inv(S)
            S_decomp = decomposition(S_pred);
            % Perform gauss elimination 
            KG = C_pred / S_decomp;
        else
            % Only use psuedo inverse when condition number much greater
            % than 1
            % Psuedo inverse is costly to compute and should not be default
            KG = C_pred * pinv(S_pred);
            % --- Alternative ----
            % Use least squares to find inverse of S_pred that solves the
            % system x * A = b
            % KG = lsqminnorm(C_pred, S_pred);
        
        end
    
        % Actualisation of predicted mean
        V_pred = (Y - mu_pred);
        X_next = m_pred + KG * V_pred;
        % Actualisation of predicted covariance
        PC = P_pred - KG * S_pred * KG.';

    else
        X_next = m_pred;
        % % ------- Constrain parameters to search space -----
        % % Constrain the parameters to [1e-8, 1]
        % X_next(3:end) = max(1e-8, min(1, X_next(3:end)));
        PC = P_pred;
        KG = [];
        V_pred = [];
    end
    % ---------- Save output for next round -----------
    output.ProcessNoiseCov = Q;
    output.MeasurementNoiseCov = R;
    output.ErrorCov = PC;
    output.dt = dt;
    output.V_pred = V_pred;
    
end

% function X_sigma = formSigmaPoints(X, P, n, lambda)
%     [N, ~] = size(X(:));
%     X_sigma = zeros(N, 2 * n+1);
%     X_sigma(:, 1) = X(:);
%     X_sigma(:, 2:n+1) = X(:) + sqrt(n + lambda) * sqrtm(P);
%     X_sigma(:, n+2:end) = X(:) - sqrt(n + lambda) * sqrtm(P);
% 
% end

function X_sigma = formSigmaPoints(X, P, n, lambda)
    [N, ~] = size(X(:));
    X_sigma = zeros(N, 2 * n + 1);
    X_sigma(:, 1) = X(:);
    % Calculate matrix square root
    if cond(P) >= 1e12
        % Matrix is ill conditioned
        P = P + eye(size(P)) * 1e-3;
    end
    L = real(sqrtm((n + lambda) * P));
    % % % % Calculate the Cholesky decomposition of the covariance matrix
    % % % % First check if positive definite
    % % % isSym = issymmetric(P);
    % % % d = eig(P);
    % % % tol = length(d) * eps(max(d));
    % % % isposdef = all(d > tol) && all(isreal(diag(P))) && isSym;
    % % % 
    % % % if isposdef
    % % %     % chol() returns the upper triangular matrix R such that R'*R = P
    % % %     % We need the lower triangular form, so we use R' or L
    % % %     L = chol((n + lambda) * P, 'lower');
    % % % else
    % % %     % Not positve definite. Make use of sqrtm...
    % % %     % Check condition number
    % % %     if cond(P) >= 1e12
    % % %         % Matrix is ill conditioned
    % % %         P = P + eye(size(P)) * 1e-3;
    % % %     end
    % % %     L = sqrtm((n + lambda) * P);
    % % % end
    % Set the non-central sigma points
    for i = 1:n
        X_sigma(:, i + 1) = X(:) + L(:, i);
        X_sigma(:, i + n + 1) = X(:) - L(:, i);
    end
end

function measuredState = measure(t, X, numStates, numMeasured)
    % Measure the first state for now...
    h = eye(numMeasured, numStates);
    measuredState = h * X;
end


