% BRUSHMODELMEX   Generate MEX-function simulateBrushModel_mex from
%  simulateBrushModel.
% 
% Script generated from project 'simulateBrushModel.prj' on 05-Mar-2025.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

clear all mex; close all; clc
% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config( "mex" );
% % cfg.DynamicMemoryAllocation = "Off";
cfg.TargetLang = "C++";
cfg.DynamicMemoryAllocationForFixedSizeArrays = true;
cfg.GenerateReport = true;
cfg.SIMDAcceleration = "Full";
cfg.OptimizeReductions = true;
cfg.EnableAutoParallelization = false;

%% Define argument types for entry-point 'simulateBrushModel'.
%%%%%%%%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_save = 1e3; % Hz
fs_sim = 1e3; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 10;
isRolling = true;

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codegen -config cfg simulateBrushModel_V2 -args model_input -nargout 2
%% Define argument types for entry-point 'simulateBrushModel'.
%%%%%%%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_save = 1e3; % Hz
fs_sim = 1e3; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 10;
isRolling = false;

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);
%%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codegen -config cfg simulateBrushModel_V2 -args model_input -nargout 2


%% Define argument types for entry-point 'shiftPressure'.
ARGS2 = cell(1,1);
ARGS2{1} = cell(10,1);
ARGS2{1}{1} = coder.typeof(uint16(0));
ARGS2{1}{2} = coder.typeof(uint16(0));
ARGS2{1}{3} = coder.typeof(single(0));
ARGS2{1}{4} = coder.typeof(0);
ARGS2{1}{5} = coder.typeof(0,[numBrushes  numBrushes]);
ARGS2{1}{6} = coder.typeof(0,[numBrushes  numBrushes]);
ARGS2{1}{7} = coder.typeof(0,[numBrushes  numBrushes]);
ARGS2{1}{8} = coder.typeof(0,[numPoints    1],[1 0]);
ARGS2{1}{9} = coder.typeof(0);
ARGS2{1}{10} = coder.typeof(0);
%% Invoke MATLAB Coder.


%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% codegen -config cfg simulateBrushModel_V2 -args ARGS3{1}

% codegen -config cfg shiftPressure -args ARGS2{1}
