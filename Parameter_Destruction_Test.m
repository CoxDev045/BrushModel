set(0, 'DefaultFigureWindowStyle', 'docked')
%% Added Noise to pressure distribution
clear; close all; clc

numBrushes  = 20;
isRolling   = true;
fs_sim      = 1e3;
fs_save     = 1e3;
t_initial   = 0;
t_final     = 40;

[model_input, sim_solution] = main(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );
v0 = model_input.v0;
omega = model_input.omega;
dt_sim = model_input.dt_sim;
Fz = model_input.Fz;
K = min(size(model_input.omega), [], 'all');

%%
animateBrushOutput(model_input, sim_solution,false, [])
%%

shift_amount_cumulative = (cumsum(v0(t_save) * dt_sim));
shift_amount = (gradient(floor(shift_amount_cumulative)) > 0);

dA = model_input.dA;

post_process = tic;
forceX              = zeros(K, model_input.LenTime_save);
forceY              = zeros(K, model_input.LenTime_save);
forceTotal          = zeros(K, model_input.LenTime_save);
forceX_medfilt      = zeros(K, model_input.LenTime_save);
forceY_medfilt      = zeros(K, model_input.LenTime_save);
forceTotal_medfilt  = zeros(K, model_input.LenTime_save);
avg_mu              = zeros(K, model_input.LenTime_save);

forceX_filt         = zeros(K, model_input.LenTime_save);
forceY_filt         = zeros(K, model_input.LenTime_save);
forceTotal_filt     = zeros(K, model_input.LenTime_save);

Y_mag2              = zeros(round(model_input.LenTime_save / 2), K);
PSD                 = zeros(round(model_input.LenTime_save / 2), K);

lgd = cell(1, K);
lgd2 = cell(1, K);

for i = 1:K
    working_data = sim_solution{i};
    P_threshold = 0.02;
    pressure_mask = working_data.PressGrid > P_threshold; %model_input.dA
    % total_valid_points = max(sum(pressure_mask( :, :), 1)); % Count valid points

    % Intergrate stress to get force
    forceX(i, :) = squeeze(sum(working_data.tauX.* dA)) ;
    forceX_medfilt(i, :) = medfilt1(forceX(i, :));
    
    forceY(i, :) = squeeze(sum(working_data.tauY.* dA)) ;
    forceY_medfilt(i, :) = medfilt1(forceY(i, :));
    
    % % %%%%%%%%%%%%%%%%%% Save forces for regression purposes %%%%%%%%%%%%%%%%%%%%
    % % X_data = cat(2, forceX(:), forceY(:));
    % % save('BrushModel_data.mat', 'v0', 'SR', 't_save', 'X_data', '-v7.3')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    forceTotal(i, :) = squeeze(sum(working_data.TotalStress.* dA)) ;
    % forceTotal = dot(working_data( :, :, 7), P_grid.save, 1) * dA;
    forceTotal_medfilt(i, :) = medfilt1(forceTotal(i, :));

    avg_mu(i, :) = squeeze(mean(working_data.mu .* pressure_mask, 1));% .* dA / 10; % Avoid division by zero
    
    isSliding = working_data.TotalStress > 0.02 * working_data.PressGrid;
end

%%
t_ind = 19e3;
plotInd = 1:1:numBrushes^2;
numElems = sqrt(max(size(plotInd)));
tauX = reshape(working_data.tauX(plotInd, t_ind), numElems, numElems);
tauY = reshape(working_data.tauY(plotInd, t_ind), numElems, numElems);
slidInd = squeeze(isSliding(plotInd, t_ind));
hasPress = squeeze(pressure_mask(plotInd, t_ind));
X = reshape(model_input.X, numElems, numElems);
Y = reshape(model_input.Y, numElems, numElems);


plotDeformationGradient(X, Y, 1e5 * tauX, zeros(size(tauY)),slidInd, hasPress, 0.5 * forceX(t_ind), 0 )
%%
colorbar_names = cell(9,1);
colorbar_names{1} = "Pressure Distribution [MPa]";
colorbar_names{2} = "Total Displacement [mm]";
colorbar_names{3} = "Total Shear Stress [MPa]";
colorbar_names{4} = "Relative X Displacement [mm]";
colorbar_names{5} = "Relative Y Displacement [mm]";
colorbar_names{6} = "\tau_x (X Shear Stress) [MPa]";
colorbar_names{7} = "\tau_y (Y Shear Stress) [MPa]";
% colorbar_names{9} = "Sliding Condition [Boolean]";
colorbar_names{8} = "Friction Coefficient";
colorbar_names{9} = 'Sliding Velocity [mm/s]';


% % Moving Average
% winSize = 10;
% b_val = (1 / winSize) * ones(1, winSize);
% a_val = 1;
% 
% for i = 1:K
%     forceTotal_filt(i, :) = filtfilt(b_val, a_val, double(forceTotal_medfilt(i, :)));
%     forceX_filt(i, :) = filtfilt(b_val, a_val, double(forceX_medfilt(i, :)));
%     forceY_filt(i, :) = filtfilt(b_val, a_val, double(forceY_medfilt(i, :)));
% end

dt_ratio = uint32(model_input.dt_ratio);

if ~isRolling
    vel = max(abs(v0(t_save)));
else
    vel = max(abs(model_input.re / 1000 * omega(t_save)));
end

for i = 1:K
    lgd{i} = sprintf("Linear Velocity of wheel v_{max} = %.1f m/s", max(model_input.re  / 1000 * omega(t_save), [], "all"));
    lgd{i+1} = sprintf("Velocity of Road element v_{max} = %.1f m/s", vel(i));
    lgd2{i} = sprintf("v_{max} = %.1f m/s", vel(i)); 
end

T = tiledlayout('horizontal');
T.Padding = "tight";
T.TileSpacing = "tight";

nexttile
hold on
plot(t_save, omega((t_save)) * model_input.re / 1000)
plot(t_save, v0((t_save)))
hold off
grid on
title('Input Velocities')
legend(lgd)
xlabel('Time[s]');ylabel('Velocity [m/s]')

nexttile
plot(t_save, forceTotal)
grid on
title('Total Force')
xlabel('Time[s]');ylabel('Force [N]')
% ylim([0, 2 * Fz])
legend(lgd2)

for i = 1:K
    [Y_mag2(:, i), half_freqz] = One_sided_fft(forceTotal(i, :), fs_save);
    N = 2 * length(Y_mag2(:, i));

    PSD(:, i) = (1 /(fs_sim * N)) * real(Y_mag2(:, i)  .* conj(Y_mag2(:, i) ));
end


figure
loglog(half_freqz, Y_mag2)

grid on
title('PSD of Force and Velocity')
xlabel('Frequency [Hz]');ylabel('Magnitude [N^2/Hz]')
legend(lgd2)

for i = 1:K
    lgd2{i} = sprintf("Simulated: v_{max} = %.1f m/s", vel(i)); 
end


figure
T = tiledlayout('horizontal');
T.Padding = "tight";
T.TileSpacing = "tight";

nexttile
plot(t_save, (forceX));
hold on
grid on;
title("Longitudinal Force [N]")
legend(lgd2)
hold off
xlabel('Time [s]');ylabel('Force [N]')
% ylim([0, 2 * Fz])

nexttile
plot(t_save, forceY);
grid on;
title("Lateral Force [N]")
legend(lgd2)
hold off
xlabel('Time [s]');ylabel('Force [N]')

for i = 1:K
    lgd2{i} = sprintf("Total Force: v_{max} = %.1f m/s",vel(i) ); 
    lgd2{i+K} = sprintf("Avg mu X F_z: v_{max} = %.1f m/s", vel(i) );
end

nexttile
plot(t_save, forceTotal);
hold on
plot(t_save, (avg_mu .* Fz), '--');
title('Total Force')
legend(lgd2)
grid on
% ylim([0, 2 * Fz])
hold off
xlabel('Time [s]');ylabel('Force [N]')

%%%%%%%%%%%%%%%%%% Plot Force vs Slip %%%%%%%%%%%%%%%%%%%%%%%
% Identify time range where slip ranges from 0% to 100%
ind = true(length(forceX), 1);%(t_save >= 11) .* (t_save <= 101);
% % Remove all the indices where ind is less than 1
% ind = ind > 0;

% Plot force vs slip using indices calculated
if isRolling
    SR_for_plot = model_input.SR(1:dt_ratio:end);
    figure
    plot(SR_for_plot(ind), forceTotal(ind))
    hold on
    plot(SR_for_plot(ind), avg_mu(ind) .* Fz)
    grid on
    xlabel('Longitudinal Slip')
    ylabel('Force Simulated [N]')
    title('Force vs Slip graph generated from Brush Model')
    legend('Simulated Force [N]', '\mu_{avg} \times F_z [N]', Location='best')
else
    v0_for_plot = model_input.v0(t_save);
    figure
    plot(v0_for_plot(ind), abs(forceX(ind)))
    hold on
    plot(v0_for_plot(ind), avg_mu(ind) .* Fz)
    grid on
    xlabel('Longitudinal Slip Velocity [mm/s]')
    ylabel('Force Simulated [N]')
    title('Force vs Slip graph generated from Brush Model')
    legend('Simulated Force [N]', '\mu_{avg} \times F_z [N]', Location='best')
end
toc(post_process)


plot_ind = 1:100:model_input.LenTime_save;
figure
T = tiledlayout('flow');
T.Padding = "tight";
T.TileSpacing = "tight";
for i = 1:K
    working_data = sim_solution{i};
    
    disp = squeeze(working_data.delta_x(1:200, plot_ind)).';
    dispY = squeeze(working_data.delta_y(1:200, plot_ind)).';
    force = squeeze(working_data.tauX(1:200, plot_ind)).';
    forceY = squeeze(working_data.tauY(1:200, plot_ind)).';
    slide_vel = squeeze(working_data.vs(1:200, plot_ind)).';
    
    % subplot(2,3,i)
    nexttile
    plot(abs(disp), force, '.')
    grid on
    % ylim([-0.2, 0.2])
    title(sprintf('Phase Plot X Stress vs X disp') )
    xlabel('Rel. Disp. [mm]')
    ylabel('Long. Stress [MPa]')

    nexttile
    plot(dispY, forceY, '.')
    grid on
    % ylim([-0.2, 0.2])
    title(sprintf('Phase Plot Y Stress vs Y disp') )
    xlabel('Rel. Disp. [mm]')
    ylabel('Lat. Stress [MPa]')
    
    
    % subplot(2,3,i+3)
    nexttile
    plot(slide_vel, force, '.')
    grid on
    % ylim([-0.2, 0.2])
    title( sprintf('Phase Plot X Stress vs X Vel' ) )
    xlabel('Sliding Vel. [mm/s]')
    ylabel('Long. Stress [MPa]')

    nexttile
    plot(slide_vel, forceY, '.')
    grid on
    % ylim([-0.2, 0.2])
    title( sprintf('Phase Plot Y Stress vs Y disp' ) )
    xlabel('Sliding Vel. [mm/s]')
    ylabel('Lat. Stress [MPa]')
end

% % figure(7)
% % comet(disp(:, 10), force(:, 10))
% % grid on
% % ylim([-0.2, 0.2])
