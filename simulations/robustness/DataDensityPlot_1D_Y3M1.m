function [ f , range90] = DataDensityPlot_1D_Y3M1( x, y, levels, fig_id, sp_id, robust_setup)
%DATADENSITYPLOT Plot the data density 
%   Makes a contour map of data density
%   x, y - data x and y coordinates
%   levels - number of contours to show
%
% By Malcolm Mclean
%   
%     if exist('limits','var')
%         map = dataDensity(x, y, 256, 256, limits);
%     else
%         map = dataDensity(x, y, 256, 256);
%        [map,range90] = dataDensity_1D(x, y, robust_setup.dpVal, robust_setup.dpVal, robust_setup);
%        [map,range90] = dataDensity_1D(x, y, robust_setup.dpVal, robust_setup.dpVal_yaxis, robust_setup);
       [map,range90] = dataDensity_2D_histogram_Y3M1(x, y, robust_setup.dpVal, robust_setup.dpVal_yaxis, robust_setup);
%     end
    map = map - min(min(map)); % takes out the 'excess' minimum value
    [num_rows, num_cols] = size(map);
% % % %     for ncol = 1:num_cols
% % % %         for nrow = 1:num_rows
% % % %             if map(nrow,ncol) ~= 0
% % % % %                 map(nrow,ncol) = map(nrow,ncol) + max(max(map))*0.05;
% % % %                 map(nrow,ncol) = map(nrow,ncol) + max(max(map))*robust_setup.shading_increase_ratio;
% % % %             end
% % % %         end
% % % %     end
% % % %     map = floor(map ./ max(max(map)) * (levels-1)); % brings into the proportion of levels, all down to zero
    for ncol = 1:num_cols
        map(:,ncol) = floor(map(:,ncol) ./ max(max(map(:,ncol))) * (levels-1)); % brings into the proportion of levels, all down to zero
    end
    map = flip(map,1); % flipped
%     figure, histogram(map)
%     figure, heatmap(map)
    
    % forcing increase for more visual
    if(isfield(robust_setup, 'forced_density')&&(robust_setup.forced_density == 1))
        % 
        [nrow,ncol] = size(map);
        for j = 1:ncol
            for i = 1:nrow
                if map(i,j) ~= 0
                    fold_factor = 1;
                    map(i,j) = map(i,j) + map(i,j) * (max(max(map)) - map(i,j)) / max(max(map)) * fold_factor;
                end
            end
        end
    end

    % set figure and position
    f = figure(fig_id);
    
    set(gca, 'Position', sp_id.Position)
    
    image(map);
% %     tempim=image(map);
%     colormap(jet(levels));
%     colormap(gray(levels));
%     colormap(flipud(gray(levels)));
    colormap(flipud(hot(levels)));
%     colormap(flipud(autumn(levels)));
    
    % setting the zero value as first in x-axis
    if robust_setup.plotSS == 1
        temp_xlim = xlim;
        set(gca, 'XLim', [0 temp_xlim(2)]);
    else
    end
    
    set(gca, 'XTick', [0 robust_setup.dpVal]);
% % %     if robust_setup.plotSS == 0
% % %         set(gca, 'XTickLabel', [min(x) max(x)]);
% % %     end
    set(gca, 'YTick', [1 robust_setup.dpVal]);
%     set(gca, 'YTickLabel', [min(y) max(y)]);
    set(gca, 'YTickLabel', [max(y) min(y)]); % flipped
%     uiwait;

    % plot range90
%     hold on
%     plot(range90(1,:), range90(3,:), 'k-')

end

