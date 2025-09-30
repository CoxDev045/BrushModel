clear;close all;clc;

Forces = [1, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6, 1, 7] * 1000;
X_len  = [98, 139, 169, 182, 200, 229, 169, 207, 225, 199, 213, 223, 101, 231];
Y_len  = [80, 141, 177, 193, 215, 234, 192, 214, 228, 182, 210, 232, 92, 235];

model = @(p) norm([p(1) * (Forces / p(2)).^(p(3)) - X_len;
                   p(4) * (Forces / p(2)).^(p(5)) - Y_len], 2);


% figure
% scatter(Forces, X_len)
% hold on
% scatter(Forces, Y_len)
% grid on
% legend('Contact Length', 'Contact Width')

K = 5;
numParticles = 12 * K;
particlePositions_unscaled = lhsdesign(numParticles, K, 'criterion','maximin');

% w_init = 10 * randn(numParticles,5) + [230, 560*9.81, 1.0, 230, 1.0];
% 
% Fmin_opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton',...
%                          'Display', 'iter',...
%                          'FiniteDifferenceType','central', ...
%                          'MaxFunctionEvaluations',3000, 'MaxIterations', 100, ...
%                          'OptimalityTolerance', 1e-10, 'StepTolerance',(eps), 'FunctionTolerance', eps, ...
%                          'UseParallel', false);
% 
% PSO_options = optimoptions('particleswarm', ...
%                            'MaxIterations', 1000, ...
%                            'SwarmSize', numParticles, ...
%                            'InitialSwarmMatrix', w_init,...
%                            'Display', 'iter', ...
%                            'UseParallel', false, ...
%                            'HybridFcn', {@fminunc, Fmin_opts}); 
% 
% % Run Optimizer
% [p_fit, gbestVal, ~, ~, ~] = particleswarm(model, K, [], [],PSO_options);

% % Initial Guess
% w_init = 2 * randn(1,5) + [230, 560*9.81, 1.0, 230, 1.0];
% 
% Fmin_opts = optimoptions('lsqnonlin', 'Algorithm', 'interior-point',...
%                          'Display', 'iter',...
%                          'FiniteDifferenceType','central', ...
%                          'MaxFunctionEvaluations',5000000, 'MaxIterations', 200000, ...
%                          'OptimalityTolerance', 1e-6, 'StepTolerance',(eps), 'FunctionTolerance', eps, ...
%                          'UseParallel', false);
% % Run Optimizer
% [p_fit, gbestVal, ~, ~, ~] = lsqnonlin(model,w_init, [], [], Fmin_opts);

% w_init = 10 * randn(numParticles,5) + [230, 560*9.81, 1.0, 230, 1.0];
% 
% OutFunc = @(options, state, flag) optimDataRecorder(options, state, flag, model);
% 
% PS_opts = optimoptions('patternsearch', ...
%                        'FunctionTolerance', eps, ...
%                        'Algorithm','nups-mads',...
%                        'Display','iter',...
%                        'MaxFunctionEvaluations',500e3,...
%                        'MaxIterations',100e3,...
%                        'SearchFcn','rbfsurrogate'); %,...
%                        % 'OutputFcn',OutFunc);
% 
% GA_options = optimoptions('ga', ...
%                            'MaxGenerations', 100, ...
%                            'PopulationSize', numParticles, ...
%                            'InitialPopulationMatrix',w_init,...
%                            'FunctionTolerance', eps,...
%                            'Display', 'iter', ...
%                            'UseParallel', false, ...
%                            'HybridFcn', {@patternsearch, PS_opts}); %,...
%                            % 'OutputFcn',OutFunc);
% % Run Optimizer
% [p_fit, gbestVal, ~, ~, ~] = ga(model, K, [], [],[], [],[], [],[],GA_options);

Fmin_opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton',...
                         'Display', 'iter',...
                         'FiniteDifferenceType','central', ...
                         'MaxFunctionEvaluations',5000, 'MaxIterations', 2000, ...
                         'OptimalityTolerance', 1e-6, 'StepTolerance',(eps), 'FunctionTolerance', eps, ...
                         'UseParallel', false);

lb = zeros(1,K);
ub = ones(1,K) * realmax;

problem = createOptimProblem('fminunc',...
                             'objective', model,...
                             'x0', [230, 560*9.81, 1.0, 230, 1.0],...
                             'lb', lb, 'ub', ub,...
                             'options', Fmin_opts);

gs = GlobalSearch('StartPointsToRun','bounds',  "PlotFcn","gsplotbestf");
% [p_fit, fval] = run(gs, problem);
ms = MultiStart(gs);
[p_fit,fval,eflag,output,manymins] = run(ms,problem,50);

% -------------------------------------------------------------
% Extract the fitted parameters
% -------------------------------------------------------------
L0 = p_fit(1);
Fz_0 = p_fit(2);
q_L = p_fit(3);
W0 = p_fit(4);
q_W = p_fit(5);

% Display the fitted parameters
fprintf('Fitted Parameters:\n');
fprintf('  L0 = %f\n', L0);
fprintf('  Fz_0 = %f\n', Fz_0);
fprintf('  q_L = %f\n', q_L);
fprintf('  W0 = %f\n', W0);
fprintf('  q_W = %f\n', q_W);

x_model = @(Fz) L0 * (Fz / Fz_0).^(q_L);
y_model = @(Fz) W0 * (Fz / Fz_0).^(q_L);

Fz_fine = linspace(100, 8000, 500);

XFine = x_model(Fz_fine);
YFine = y_model(Fz_fine);

figure

nexttile
plot(Fz_fine, XFine)
hold on
scatter(Forces, X_len)
ylabel('Length [m]')
xlabel('Vertical Load [N]')

nexttile
plot(Fz_fine, YFine)
hold on
scatter(Forces, Y_len)
ylabel('Length [m]')
xlabel('Vertical Load [N]')

