set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear; close all; clc

numBrushes  = floor(linspace(4, 90, 15 ));
isRolling   = true;
fs_sim      = 1e3;
fs_save     = 1e3;
t_initial   = 0;
t_final     = 40;
time_to_sim = zeros(size(numBrushes));

for i = 1:length(numBrushes)
    fprintf('Starting sim in a grid of [%d x %d] brushes\n', numBrushes(i), numBrushes(i))
    start = tic;
    [model_input{i}, sim_solution{i}] = main(numBrushes(i), isRolling, fs_sim, fs_save, t_initial, t_final);
    time_to_sim(i) = toc(start);
end

K = length(numBrushes);
% whos sim_solution omega v0
%%
save("Mesh_dependency.mat", "model_input", "sim_solution", "numBrushes", "time_to_sim", '-v7.3')
%%
post_process = tic;
forceX              = zeros(K, model_input{1}.LenTime_save);
forceY              = zeros(K, model_input{1}.LenTime_save);
forceTotal          = zeros(K, model_input{1}.LenTime_save);

avg_mu              = zeros(K, model_input{1}.LenTime_save);

lgd = cell(1, K);
lgd2 = cell(1, K);

t_save = single( linspace(t_initial, t_final, t_final * fs_save + 1) );

for i = 1:K
    dA = model_input{1, i}.dA;
    v0 = model_input{1, i}.v0;
    omega = model_input{1, i}.omega;
    dt_sim = model_input{1, i}.dt_sim;
    Fz = model_input{1, i}.Fz;
    
    shift_amount_cumulative(i, :) = (cumsum(v0(t_save) * dt_sim));
    shift_amount(i, :) = (gradient(floor(shift_amount_cumulative(i, :))) > 0);
   
    P_threshold = 0.02;
    pressure_mask = sim_solution{1, i}{1,1}.PressGrid > P_threshold; %model_input.dA
    
    % Intergrate stress to get force
    forceX(i, :) = squeeze(sum(sim_solution{1, i}{1,1}.tauX.* dA)) ;
    forceY(i, :) = squeeze(sum(sim_solution{1, i}{1,1}.tauY.* dA)) ;
    
    % % %%%%%%%%%%%%%%%%%% Save forces for regression purposes %%%%%%%%%%%%%%%%%%%%
    % % X_data = cat(2, forceX(:), forceY(:));
    % % save('BrushModel_data.mat', 'v0', 'SR', 't_save', 'X_data', '-v7.3')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    forceTotal(i, :) = squeeze(sum(sim_solution{1, i}{1,1}.TotalStress.* dA)) ;

    avg_mu(i, :) = squeeze(mean(sim_solution{1, i}{1,1}.mu .* pressure_mask, 1));

    lgd{i} = sprintf('Grid size: [%d x %d]', numBrushes(i), numBrushes(i));
end
%%
figure
plot(t_save, forceTotal)
xlabel('Time [s]')
ylabel('Magnitude [N]')
legend(lgd)
title("Mesh dependency Study: Total Force ")
grid on

figure
plot(t_save, avg_mu .* Fz)
xlabel('Time [s]')
ylabel('Magnitude [N]')
legend(lgd)
title("Mesh dependency Study: Force from Friction Coeff. ")
grid on

figure
plot(numBrushes, time_to_sim)
xlabel('Grid size')
ylabel('Time to solve [s]')
title('Grid size vs Time to solve Mesh dependency imulation')
grid on


%%
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