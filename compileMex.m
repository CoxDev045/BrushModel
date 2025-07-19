function compileMex(model_input)
    % Create configuration object of class 'coder.MexCodeConfig'.
    cfg = coder.config( "mex" );
    % % cfg.DynamicMemoryAllocation = "Off";
    cfg.TargetLang = "C++";
    cfg.DynamicMemoryAllocationForFixedSizeArrays = true;
    cfg.GenerateReport = true;
    cfg.SIMDAcceleration = "Full";
    cfg.OptimizeReductions = true;
    cfg.EnableAutoParallelization = false;
        
    %%%%%%%%%% Rolling Tyre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Invoking MATLAB Coder now...\n')
    try
        codegen -config cfg simulateBrushModel_V2 -args model_input -nargout 2
    catch ME
        fprint('Code generation failed...\n')
        rethrow(ME)
    end

    fprintf('Successfully generated code!\n')
    S = load('gong.mat');
    sound(S.y, 2 * S.Fs)

end

