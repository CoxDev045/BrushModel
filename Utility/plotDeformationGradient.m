function plotDeformationGradient(X, Y, tauX, tauY, is_sliding, pressure_mask, Fx, Fy)
%PLOTDEFORMATIONGRADIENT Summary of this function goes here
%   Detailed explanation goes here

% Generate some sample stress components (replace with your sigma_x, sigma_y)
% For demonstration, let's make a simple flow pattern
sigma_x = tauX; % Stress in x depends on y
sigma_y = tauY;  % Stress in y depends on x

% --- 2. Create the Quiver Plot ---
figure;
hold on;

% Plot the stress vectors using quiver
% The last argument (0.8 here) scales the arrows. Adjust for better visualization.
quiver(X, Y, sigma_x, sigma_y, 0.8, 'Color', 'blue');

% Simple rectangular regions for illustration:
% Green region (assuming your simulation area is roughly a rectangle for now)
x_min = min(X(:)); x_max = max(X(:));
y_min = min(Y(:)); y_max = max(Y(:));

% Find the x and y coordinates where the brushes are sliding
isActiveAndSliding = logical(is_sliding .* pressure_mask);
x_sliding = X(isActiveAndSliding);
y_sliding = Y(isActiveAndSliding);

% Green region (above threshold)
% --- Plot the Sliding Brushes as a Green Scatter Plot ---
% Plot the coordinates of the sliding brushes with a large, green marker
scatter(x_sliding, y_sliding, 100, 'green', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5); 


% Red region (below threshold)
% --- Plot the Sticking Brushes as a Red Scatter Plot ---
% You can also plot the sticking brushes with a different color
isActiveAndSticking = logical(~is_sliding .* pressure_mask);
x_sticking = X(isActiveAndSticking);
y_sticking = Y(isActiveAndSticking);
scatter(x_sticking, y_sticking, 100, 'red', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);

% Red region (below threshold)
% --- Plot the Sticking Brushes as a Red Scatter Plot ---
% You can also plot the sticking brushes with a different color
isNotActive = ~pressure_mask;
x_InActive = X(isNotActive);
y_InActive = Y(isNotActive);
scatter(x_InActive, y_InActive, 100, [211,211,211] / 255, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1.0);


% --- 4. Customizing the Plot ---
axis equal; % Important to prevent distortion of arrows
xlabel('X-coordinate');
ylabel('Y-coordinate');
title('Stress Distribution');
grid on;
box on;
xlim([x_min-0.1 x_max+0.1]);
ylim([y_min-0.1 y_max+0.1]);

% You can add the specific purple arrow separately if it represents
% a resultant force or a specific point stress.
% For example, a single, larger arrow at the center:
quiver(0, 0, Fx, Fy, 0.8, 'Color', 'red', 'LineWidth', 2, 'MaxHeadSize', 0.5);
plot(0, 0, 'k.', 'MarkerSize', 15); % Black dot at the center

hold off;
legend('Tangential Stress Vectors', 'Sliding region', 'Adhesion region', 'Inactive Brushes', 'Resultant Force [kN]')
end

