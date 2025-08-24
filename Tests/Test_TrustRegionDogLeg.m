clear;close all;clc

R = 0:.002:1;
TH = 2*pi*(0:.002:1); 
X = R'*cos(TH);
Y = R'*sin(TH);
t = 1;

Z = log(1 + rosenbrock(t, [X(:),Y(:)]));
Z = reshape(Z,size(X));

figure
surf(X, Y, Z)
shading flat

max_iters = 100;
options.JacTol = 1e-3;
options.FARGS = {};
options.ThreshScal = [1;1];

x0 = rand(1, 2);

for i = 1:max_iters
    [x_next, Fty_next] = solveTrustRegionDogLeg(@rosenbrock, t, x0, options);
end

%%
function y = rosenbrock(t, x, varargin)
[~, M] = size(x);
if M ~= 2
    % Input is column vector
    xp = x(:).';
else
    xp = x;
end

y = 100 * (xp(:, 2) - xp(:, 1).^2).^2 + (1 - xp(:, 1)).^2;
end