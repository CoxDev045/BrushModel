function J = approximateJacobianComplex(t, w, func, args)
    % Compute numerical Jacobian of a function using complex-step differentiation
    % J has size D x K for each time step (D states, K parameters)
    
    % Input size
    K = numel(w); 
    
    % Get function output and its size
    x = func(t, w, args{:});
    D = numel(x); 
    
    % Initialize Jacobian with the correct dimensions: (output size) x (input size)
    J = zeros(D, K);
    
    h = 1e-20; % A very small number for complex step
    
    % Loop through each input parameter to perturb it
    for i = 1:K
        w_complex = w;
        
        % Perturb the i-th input variable with a complex step
        w_complex(i) = w_complex(i) + 1i * h;
        
        % Evaluate the function with the complex input
        x_complex = func(t, w_complex, args{:});
        
        % The imaginary part of the result, divided by h, gives the i-th column of the Jacobian
        J(:, i) = imag(x_complex) / h;
    end
end