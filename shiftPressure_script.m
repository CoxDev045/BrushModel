% SHIFTPRESSURE_SCRIPT   Generate MEX-function shiftPressure_mex from
%  shiftPressure.
% 
% Script generated from project 'shiftPressure.prj' on 10-Mar-2025.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.EnableRuntimeRecursion = false;
cfg.GenerateReport = true;
cfg.NumberOfCpuThreads = 6;
cfg.ReportPotentialDifferences = false;
cfg.SIMDAcceleration = 'Full';
cfg.EnableAutoParallelization = true;
cfg.OptimizeReductions = true;

%% Define argument types for entry-point 'shiftPressure'.
ARGS = cell(1,1);
ARGS{1} = cell(10,1);
ARGS{1}{1} = coder.typeof(uint16(0));
ARGS{1}{2} = coder.typeof(uint16(0));
ARGS{1}{3} = coder.typeof(single(0));
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0,[400  1]);
ARGS{1}{6} = coder.typeof(0,[400  1]);
ARGS{1}{7} = coder.typeof(0,[400  1]);
ARGS{1}{8} = coder.typeof(0,[1000000    1],[1 0]);
ARGS{1}{9} = coder.typeof(0);
ARGS{1}{10} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg shiftPressure -args ARGS{1}

