% BRUSHMODELMEX   Generate MEX-function simulateBrushModel_mex from
%  simulateBrushModel.
% 
% Script generated from project 'simulateBrushModel.prj' on 05-Mar-2025.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

% SHIFTPRESSURE_SCRIPT   Generate MEX-function shiftPressure_mex from
%  shiftPressure.
% 
% Script generated from project 'shiftPressure.prj' on 07-Mar-2025.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

clear; close all; clc
%% Create configuration object of class 'coder.MexCodeConfig'.
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
numBrushes = 20;
numNodes = numBrushes^2;
fs_sim = 1e5;
time_final = 10;
numPoints = fs_sim * time_final + 1;

ARGS = cell(1,1);
ARGS{1} = cell(13,1);
ARGS{1}{1} = coder.typeof(uint16(0));
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(single(0),[numNodes numPoints]);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0,[numPoints    1]);
ARGS{1}{8} = coder.typeof(0);
ARGS{1}{9} = coder.typeof(0);
ARGS{1}{10} = coder.typeof(0);
ARGS{1}{11} = coder.typeof(0);
ARGS{1}{12} = coder.typeof(0,[numNodes  1]);
ARGS{1}{13} = coder.typeof(0,[numNodes  1]);

%% Define argument types for entry-point 'simulateBrushModel'.
%%%%%%%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numBrushes = 20;
numNodes = numBrushes^2;
fs_sim = 1e5;
time_final = 3;
numPoints = fs_sim * time_final + 1;

ARGS3 = cell(1,1);
ARGS3{1} = cell(12,1);
ARGS3{1}{1} = coder.typeof(uint16(0));
ARGS3{1}{2} = coder.typeof(0);
ARGS3{1}{3} = coder.typeof(0);
ARGS3{1}{4} = coder.typeof(single(0),[numNodes numPoints]);
ARGS3{1}{5} = coder.typeof(0);
ARGS3{1}{6} = coder.typeof(0);
ARGS3{1}{7} = coder.typeof(0,[numPoints    1]);
ARGS3{1}{8} = coder.typeof(0);
ARGS3{1}{9} = coder.typeof(0);
ARGS3{1}{10} = coder.typeof(0);
ARGS3{1}{11} = coder.typeof(0,[numNodes  1]);
ARGS3{1}{12} = coder.typeof(0,[numNodes  1]);


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
%%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% codegen -config cfg simulateBrushModel_V2 -args ARGS{1}

%%%%%%%%% Sliding Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codegen -config cfg simulateBrushModel_V2 -args ARGS3{1}

% codegen -config cfg shiftPressure -args ARGS2{1}
