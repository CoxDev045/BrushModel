function [state,options,optchanged] = optimDataRecorder(options, state, flag, objective)
    persistent h1 history s
    optchanged = false;
    switch flag
        case 'init'
            h1 = figure;
            T = tiledlayout('flow');
            T.Padding = "compact";
            T.TileSpacing = 'tight';
            s = scatter(state.Generation,  log10(min(state.Score)));
            grid on
            xlabel('Generation')
            ylabel('Score')
            history(:,:,1) = state.Population;
            assignin('base','gapopulationhistory',history);
        case 'iter'
            % Update the history every 10 generations.
            if rem(state.Generation,10) == 0
                ss = size(history,3);
                history(:,:,ss+1) = state.Population;
                assignin('base','gapopulationhistory',history);
            end
            % Find the best objective function, and stop if it is low.
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);
            bestf = objective(bestx);
            if bestf <= 0.1
                state.StopFlag = 'y';
                disp('Got below 0.1')
            end
            % Update the plot.
            figure(h1)
            s.XData = [s.XData, state.Generation]; s.YData = [s.YData, log10(min(state.Score))];
            % pause(0.001)
            % Update the fraction of mutation and crossover after 25 generations.
            if state.Generation == 25
                options.CrossoverFraction = 0.8;
                optchanged = true;
            end
        case 'done'
            % Include the final population in the history.
            ss = size(history,3);
            history(:,:,ss+1) = state.Population;
            assignin('base','gapopulationhistory',history);
    end
end