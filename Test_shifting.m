clear; close all; clc

[X, Y] = meshgrid(linspace(-100, 20, 20), linspace(-100,20, 20));
Press = 0.2 * sin((X + Y) * 0.2) + 0.2;

% surf(X, Y, Press)
PressX = X; PressY = Y;
X = X(:); Y = Y(:);%Press = Press(:);


maxX = max(X, [], 'all');
maxY = max(Y, [], 'all');
minX = min(X, [], 'all');
minY = min(Y, [], 'all');

wrap_range = maxX - minX;

tempX = mod(X, wrap_range) + minX;

shift_amount = cumsum(linspace(0, 10, 100));

X_shifted = zeros(length(X), length(shift_amount));
tempPress = zeros(length(PressX), length(PressY), length(shift_amount));

for i = 1:length(shift_amount)
    X_shifted(:, i) = mod(tempX + shift_amount(i), wrap_range) + minX;

    tempPress(:, :, i) = interp2(PressX, PressY, Press, reshape(X_shifted(:, i), 20, 20), PressY, 'cubic');
end


% figure
% % subplot(211)
% plot(X)
% hold on
% plot(X_shifted, '-o')
% grid on
% xlabel('Coordinate Index')
% ylabel('Coordinate value')
% legend(arrayfun(@(i) sprintf('Shift = %.2f', shift_amount(i)), 1:length(shift_amount), 'UniformOutput', false) )

% subplot(212)
% plot(shift_amount)

%

% Plot an initial frame and store the handle to the image object
im = surf(PressX, PressY, tempPress(:, :, 1));  % Initial frame
shading interp
c = colorbar;
% c.Label.String = "X values [mm]";
title("Pressure Distribution over time");
ylabel('Lateral y-direction [mm]');
xlabel('Longitudinal x-direction [mm]');

zlim([0, 0.4])
clim([0, 0.4])


plot_ind = 1:length(shift_amount);

% Loop to update the image data
for t = plot_ind
    % Update the data in the existing image plot rather than redrawing it
    set(im, 'ZData', tempPress(:, :, t));  % Update the data
    
    pause(0.1)
end
