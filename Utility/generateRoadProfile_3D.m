function [roadProfile] = generateRoadProfile_3D(class, spatial_sampling_frequency_x,spatial_sampling_frequency_y, roadLen, roadWidth)
    % % Function to generate the 2D spatial representation of a specified road class. 
    % % Method was extended from the on defined in the ISO 8608:1995 standard.
    % % Authors: Devon Cox; Andries J. Peenze 
    % % Date: July 2025 
    % %
    % % class                         % Class of road as defined in ISO 8608:1995
    % %
    % % spatial_sampling_frequency_x  % Sampling frequency in the X-direction of the grid (cycles / metre)
    % % spatial_sampling_frequency_y  % Sampling frequency in the Y-direction  of the grid (cycles / metre)
    % % 
    % % roadLen     % Distance of road in x direction (metres)
    % % roadWidth   % Width of road in y direction (metres)
    
    num_points_x = roadLen * spatial_sampling_frequency_x; % Number of points in X-direction
    num_points_y = roadWidth * spatial_sampling_frequency_y; % Number of points in Y-direction

    freq_res_x = spatial_sampling_frequency_x / num_points_x;
    freq_res_y = spatial_sampling_frequency_y / num_points_y;
    
    
    % --- 1. Create the X and Y coordinate grids ---
    % Create 1D coordinate vectors
    x_coords = linspace((1e-4), (spatial_sampling_frequency_x), num_points_x);
    y_coords = linspace((1e-4), (spatial_sampling_frequency_y), num_points_y);
    
    
    % Create 2D meshgrid for coordinates
    % X_grid: rows are constant X, columns vary
    % Y_grid: rows vary Y, columns constant
    [X_grid, Y_grid] = meshgrid(x_coords, y_coords);
    
    % --- 2. Determine the center of the grid ---
    x_center = spatial_sampling_frequency_x / 2;
    y_center = spatial_sampling_frequency_y / 2;
    
    % --- 3. Calculate the radial distance from the center for each point ---
    % cone_Radius_x = spatial_sampling_frequency_x / 2;
    % cone_Radius_y = spatial_sampling_frequency_y / 2;
    % Maximum value that R_grid may be
    Max_allowed_radius = x_center;%hypot(cone_Radius_x, cone_Radius_y);
    Min_allowed_radius = 1;
        
    % Evaluate x and y grid on radial function
    n = 2;
    R_grid = radialFunction(X_grid, x_center, ...
                            Y_grid, y_center, ...
                            n, Max_allowed_radius, Min_allowed_radius);



    
    % Get the half indices for further calcs
    [N, M] = size(R_grid);
    halfIndX = calculateInd(N);
    halfIndY = calculateInd(M);
    
    %%%%%%%%%%% Method found in literature %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Reference: Kropáč, O. and Múčka, P., 2009. 
    % Classification scheme for random longitudinal road unevenness considering road waviness and vehicle response.
    % Shock and Vibration, 16(3), pp.273-289.
    %
    
    if ischar(class)
        availableClasses = 'ABCDEFGH';
        selectedClass = find(availableClasses(:) == class);
        
        n0 = 1/(2*pi);
        alpha = 1 + selectedClass;
        Gd_n0 = (4)^(alpha) * 1e-6;
        w = 2;
        % Generate PSD profile of 1D road profile evaluated on distance
        % function grid
        z = psdRoadProfile(n0, Gd_n0, w, R_grid);
    elseif isnumeric(class)
        w = 2; % Assume a gradient of -2
        MaxVal = class;
        % Generate PSD Profile based on MaxVal supplied
        z = psdSurfaceProfile(MaxVal, w, R_grid);
    else
        error('Please supply a character for road input or a maximum value for surface input!')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For IFFT2 the spikes have to be at the four corners and not at the middle
    % of the plane
    topLeft = z(    halfIndX:-1:1,     halfIndY:-1:1);
    botLeft = z(    end:-1:halfIndX+1,   halfIndY:-1:1);
    topRight = z(   halfIndX:-1:1,     end:-1:halfIndY+1);
    botRight = z(   end:-1:halfIndX+1,   end:-1:halfIndY+1);
    
    tempZ = [topLeft, topRight;
             botLeft, botRight];
    
    % Scale PSD values to reflect those of what FFT would have been
    fftMag = sqrt( tempZ * 2 * ( freq_res_x * freq_res_y) ) * ( N * M);
    
    phaseAngle = rand(N, M) * 2 * pi - pi;

    % Ensure DC values (Corners) has zero phase
    phaseAngle(1, 1) = 0;
    phaseAngle(end, 1) = 0;
    phaseAngle(1, end) = 0;
    phaseAngle(end, end) = 0;
    
    % Ensure Nyquist values (middle of edges) has zero phase
    phaseAngle(1, halfIndY) = 0; % Middle of first row
    phaseAngle(halfIndX, 1) = 0; % Middle of first column
    phaseAngle(halfIndX, halfIndY) = 0; % Middle of grid
    phaseAngle(end, halfIndY) = 0; % Middle of last row
    phaseAngle(halfIndX, end) = 0; % Middle of last column

    % Apply phase angles to magnitude
    z = fftMag .* exp(1j * phaseAngle);
    
    % Ensure DC are purely real
    z(1, 1) = real(z(1, 1));
    z(end, 1) = real(z(end, 1));
    z(1, end) = real(z(1, end));
    z(end, end) = real(z(end, end));
   
    % Ensure Nyquist values (middle of edges) are purely real
    z(1, halfIndY)         = real(z(1, halfIndY)       ); % Middle of first row
    z(halfIndX, 1)         = real(z(halfIndX, 1)       ); % Middle of first column
    z(halfIndX, halfIndY)  = real(z(halfIndX, halfIndY)); % Middle of grid
    z(end, halfIndY)       = real(z(end, halfIndY)     ); % Middle of last row
    z(halfIndX, end)       = real(z(halfIndX, end)     ); % Middle of last column


    % Take inverse FFT to translate back to spatial domain. The symmetric
    % flag ensures the output is real. See documentation for further
    % details
    roadProfile = ifft2(z, N, M, 'symmetric') / 1e2;

end

function z = radialFunction(x, xc, y, yc, n, Max_allowed_radius, Min_allowed_radius)

    z = ((x - xc).^n + (y - yc).^n).^(1/n);
    z = min(z, Max_allowed_radius);
    z = max(z, Min_allowed_radius);
    
end

function halfInd = calculateInd(N)

    if mod(N, 2) == 0
        halfInd = round(N / 2);
    else
        halfInd = round(N/2) + 1;
    end
end

function y = psdRoadProfile(n0, Y0, w, n)

    y = Y0 * (n / n0).^(-w);
end

function y = psdSurfaceProfile(MaxVal, w, n)

    y = MaxVal * (n).^(-w);
end


