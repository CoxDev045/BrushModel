clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear; close all; clc

folder = 'Slip_dependency_study'; % Path to saved results
files = dir(fullfile(folder, 'Slip_dependency_study_n*_SR*.mat'));

% Preallocate storage for results
numFiles = length(files);
results(numFiles) = struct('n', [], 'SR', [], 'time_span', [], 'forceX', [], 'forceY', [], 'forceTotal', [], 'avg_mu', []);

validCount = 0;
for i = 1:numFiles
    filePath = fullfile(folder, files(i).name);
    
    % Extract n (mesh size) and SR from the filename
    tokens = regexp(files(i).name, 'n(\d+)_SR([-+]?[0-9]*\.?[0-9]+)', 'tokens');
    if isempty(tokens), continue; end  % Skip if no match found
    
    validCount = validCount + 1;
    n_val = str2double(tokens{1}{1});
    SR_val = str2double(tokens{1}{2});
    
    % Load the data
    data = load(filePath);
    
    % Store data in structure
    results(validCount).n = n_val;
    results(validCount).SR = SR_val;
    results(validCount).time_span = data.time_span;
    results(validCount).forceX = data.forceX;
    results(validCount).forceY = data.forceY;
    results(validCount).forceTotal = data.forceTotal;
    results(validCount).avg_mu = data.avg_mu;
end

% Remove empty entries
results = results(1:validCount);

[~, sortIdx] = sortrows([[results.n]', [results.SR]'], [1, 2]); % Sort by n, then SR
results = results(sortIdx);

unique_n = unique([results.n]);
time_len = floor(length(results(1).time_span)/2);

figure; hold on;
for j = 1:length(unique_n)
    idx = [results.n] == unique_n(j);
    tempForce = {results(idx).forceTotal};
    meanForces = cellfun(@(x) mean(x(time_len:end)), tempForce, 'UniformOutput', true);

    SR_values = [results(idx).SR];  % Extract SR values as an array

    meanForces = arrayfun(@(sr, force) force * (-1)^(sr < 0), SR_values, meanForces);
    
    plot(SR_values, meanForces, '-o', 'DisplayName', sprintf('n = %d', unique_n(j)^2));
end
xlabel('Slip Ratio (SR)');
ylabel('Total Force');
title('Force vs. SR for Different Mesh Sizes');
grid on;

% Plot Savkoor Friction Model

SR = -0.3:0.05:0.3;
vs = linspace(-50, 50, numel(SR));
mu_s = Savkoor_Friction(vs);
mu_s(1:floor(numel(vs) / 2)) = -mu_s(1:floor(numel(vs) / 2));


figure(1)
plot(SR, 25 * mu_s, 'k--', 'DisplayName','Expected')
legend('Location','best');
% legend(lgd)
% hold on
% plot([23.47], [1.15], 'ro', 'MarkerFaceColor','r')
% plot(-[23.47], -[1.15], 'ro', 'MarkerFaceColor','r')
% grid on
% xlabel('Sliding Velocity [m/s]')
% ylabel('Friction COefficient')
% legend('', 'Max Friction Coefficient')

%%
figure; hold on;
for j = 1:length(unique_n)
    idx = [results.n] == unique_n(j);
    plot([results(idx).SR], cellfun(@mean, {results(idx).avg_mu}), '-o', 'DisplayName', sprintf('n = %d', unique_n(j)));
end
xlabel('Slip Ratio (SR)');
ylabel('Average Friction Coefficient (\mu)');
title('μ vs. SR for Different Mesh Sizes');
legend;
grid on;

%%
[N_mesh, SR_mesh] = meshgrid(unique([results.n]), unique([results.SR]));
Force_mesh = NaN(size(N_mesh));

for i = 1:numel(N_mesh)
    idx = ([results.n] == N_mesh(i) & [results.SR] == SR_mesh(i));
    if any(idx)
        Force_mesh(i) = mean(results(idx).forceTotal);
    end
end

figure;
surf(N_mesh, SR_mesh, Force_mesh);
xlabel('Mesh Size (n)');
ylabel('Slip Ratio (SR)');
zlabel('Total Force');
title('Total Force vs. Mesh Size and SR');
colorbar;
grid on;

%%
lgd = cell(numel(n)+1,1);
lgd{end} = 'Expected';

figure(1)
for i = 1:numel(n)
    plot(SR, force)
    lgd{i} = strcat(num2str(n(i)^2), ' Elements');
end
hold on
% plot([SR(1), SR(end)], -[0.02*100, 0.02*100], 'k--')
% hold off
grid on
xlabel('Slip [%]')
ylabel('Longitudinal Force')
title('Mesh Dependency Study on Force vs Slip Graph')


%%

function mu_s = Savkoor_Friction(vs)
mu_0 = 0.02;
mu_m = 1.15;
h = 0.23;
v_m = 23.47;
mu_s = mu_0 + (mu_m - mu_0) * exp(-(h * log10(abs(vs) / v_m)).^2);
% if vs < 0
%     mu_s = mu_0 + (mu_m - mu_0) * exp(-(h * log10(vs / v_m) ).^2);
% else
%     mu_s = mu_0 + (mu_m - mu_0) * exp(-(h * log10(abs(vs) / v_m)).^2);
%     mu_s = -1 * mu_s;
% end

end