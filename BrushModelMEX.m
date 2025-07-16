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

% Define argument types for entry-point 'simulateBrushModel'.
%%%%%%%%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_save = 1e3; % Hz
fs_sim = 1e4; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 120;
isRolling = true;

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codegen -config cfg simulateBrushModel_V2 -args model_input -nargout 2

load gong.mat
sound(y, 5 * Fs)

run('PostProcess_Simdata.m')

%%
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

% Define argument types for entry-point 'simulateBrushModel'.
%%%%%%%%%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_save = 1e3; % Hz
fs_sim = 1e4; % Hz
numBrushes = 20;
t_initial = 0;
t_final = 120;
isRolling = false;

model_input = brush_init(numBrushes, isRolling, fs_sim, fs_save, t_initial, t_final);

%%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codegen -config cfg simulateBrushModel_V2 -args model_input -nargout 2

load gong.mat
sound(y, 5 * Fs)

run('PostProcess_Simdata.m')
