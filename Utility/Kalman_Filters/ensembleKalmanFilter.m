function [X_next, KG, P, output] = ensembleKalmanFilter(initial_state_mean, initial_covariance, num_ensemble_members, Q, R, measurements)
    % ENSEMBLEKALMANFILTER implements the EnKF for state estimation
    %   This function simulates the EnKF over a series of measurements.
    %
    %   Inputs:
    %   initial_state_mean:     Initial guess for the state vector mean.
    %   initial_covariance:     Initial state covariance matrix.
    %   num_ensemble_members:   Number of members in the ensemble.
    %   Q:                      Process noise covariance matrix.
    %   R:                      Measurement noise covariance matrix.
    %   measurements:           A matrix of measurements over time.
    %
    %   Outputs:
    %   X_next:              The estimated state mean at each time step.
    %   KG:                  The Kalman Gain matrix at each time step.
    %   P:                   The estimated covariance at each time step.
    %   output:              
    
    % -------------------- Main Filter Loop --------------------
        % Get the current measurement
        y_k = measurements(:, k);
        
        % ================ Prediction Step =================
        % Propagate each ensemble member through the non-linear dynamics
        % function. Add process noise to each member.
        
        
        % Calculate the predicted state mean from the new ensemble
        
        
        % ================ Update Step =================
        % Calculate the predicted measurement for each ensemble member.
        % Add measurement noise to the actual measurement for the stochastic update.
        
        
        % Calculate the predicted measurement covariance from the ensemble.
        
        
        % Calculate the cross-covariance between the state and the measurements.
        
        
        % Calculate the Kalman gain.
        
        
        % Apply the stochastic update equation to each ensemble member.
        
        
        % ================ Final Step =================
        % Calculate the final state mean (the filter's output) from the updated ensemble.
        
        
        % Calculate the final covariance matrix from the updated ensemble.
        
        
        % Store the results for this time step.
        

    
end