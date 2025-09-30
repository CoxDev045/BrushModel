function animateBrushOutput(model_input, sim_solution,saveVideo, VideoOpts)
%ANIMATEBRUSHOUTPUT Summary of this function goes here
    plot_solution = struct();
    
    working_data = sim_solution{1};
    savedFields = fieldnames(working_data);

    numBrushes = model_input.numElems;
    
    for i = 1:length(savedFields)
        plot_solution.(savedFields{i}) = reshape(working_data.(savedFields{i}), numBrushes, numBrushes, model_input.LenTime_save);
    end

    fprintf('Simulation Solution is ready for animation! \n');
        
    if saveVideo
        VideoName = VideoOpts.VideoName;
        VideoDir = VideoOpts.Dir;
        FrameRate = VideoOpts.FrameRate;
        % Initialize Video
        % video_filename = sprintf('NoisyPress_10s_%dN_BrushSim_%dHz_%drpm_slip%.2f_omegaZ%.2f_alpha%.2f.mp4', Fz, fs_sim, rpm, SR, omega_z, alpha);
        video_filename = sprintf(VideoName);
        Video_path = fullfile(strcat(VideoDir,video_filename));
        v = VideoWriter(Video_path, 'MPEG-4');
        v.FrameRate = FrameRate;  
        open(v);
    end
    
    figure(7);
    h = gobjects(1, 9);
    im = gobjects(1, 9);
    
    clim_values = [
        0, 0.4;
        0, 100.0;
        0, 1.0;
        -100.0, 100.0;
        -100.0, 100.0;
        -1.0, 1.0;
        -0.5, 0.5;
        0, 3.5;
        0, 30.0;
    ];
    
    % Initialize subplots
    for j = 1:length(savedFields)
        h(j) = subplot(3, 3, j);
        im(j) = imagesc(model_input.X(1, :), model_input.Y(:, 1), plot_solution.(savedFields{j})(:, :, 1));
        c = colorbar;
        c.Label.String = savedFields{j};
        title(savedFields{j});
        ylabel('Lateral y-direction [mm]');
        xlabel('Longitudinal x-direction [mm]');
        set(h(j), 'CLim', clim_values(j, :));
    end
    
    pause(1);
    
    plot_ind = 1:50:model_input.LenTime_save;
    
    % Animation Loop
    for t = plot_ind    
        for j = 1:9
            set(im(j), 'CData', plot_solution.(savedFields{j})(:, :, t));
        end
    
        % Capture frame
        frame = getframe(gcf);
        % % writeVideo(v, frame);
    
        % Reduce lag
        if mod(t, 10) == 0
            pause(0.001);
        end
    end
    
    % Finalize Video
    if saveVideo
        close(v);
    end
    fprintf("Animation Successfully saved as Video!\n");

end

