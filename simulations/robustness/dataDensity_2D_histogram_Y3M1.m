% % % % function [ dmap ] = dataDensity_1D( x, y, width, height, limits, fudge )
function [ dmap,range90 ] = dataDensity_2D_histogram_Y3M1( x, y, width, height, robust_setup)
%DATADENSITY Get a data density image of data 
%   x, y - two vectors of equal length giving scatterplot x, y co-ords
%   width, height - dimensions of the data density plot, in pixels
%   limits - [xmin xmax ymin ymax] - defaults to data max/min
%   fudge - the amount of smear, defaults to size of pixel diagonal
%

% % % % % Safecopy initial code
% % % % % By Malcolm McLean
% % % % %
% % % % %     if(nargin == 4)
% % % %         limits(1) = min(x);
% % % %         limits(2) = max(x);
% % % %         limits(3) = min(y);
% % % %         limits(4) = max(y);
% % % % %     end
% % % %     deltax = (limits(2) - limits(1)) / width;
% % % %     deltay = (limits(4) - limits(3)) / height;
% % % % %     deltay = (limits(3) - limits(4)) / height;
% % % %     if(nargin < 6)
% % % %         fudge = sqrt(deltax^2 + deltay^2);% / 1;
% % % % %         fudge = 1;
% % % %     end
% % % %     dmap = zeros(height, width);
% % % %     for ii = 0: height - 1
% % % %         yi = limits(3) + ii * deltay + deltay/2;
% % % %         for jj = 0 : width - 1
% % % %             xi = limits(1) + jj * deltax + deltax/2;
% % % %             dd = 0;
% % % %             for kk = 1: length(x)
% % % %                 dist2 = (x(kk) - xi)^2 + (y(kk) - yi)^2; % here the distance to each point in the loop is calculated.
% % % % %                 dist2 = (x(kk) - xi)^2; % here the distance to each point in the loop is calculated.
% % % % %                 dist2 = (x(kk) - xi)^2.05 + (y(kk) - yi)^2.05; % here the distance to each point in the loop is calculated.
% % % % %                 dist2 = 4*(x(kk) - xi)^2 + (y(kk) - yi)^2; % here the distance to each point in the loop is calculated.
% % % % %                 dist2 = (x(kk) - xi)^3 + (y(kk) - yi)^3; % here the distance to each point in the loop is calculated.
% % % %                 dd = dd + 1 / ( dist2 + fudge); 
% % % %             end
% % % %             dmap(ii+1,jj+1) = dd;
% % % %         end
% % % %     end
% % % %     range90 = 0;
            

% % % %     % Update version. Time series, 1-dimensional.
% % % %     % 1a. Create empty 90% range array
% % % %     dmap = zeros(height+1, width+1);
% % % %     % 1b. Create array with complete size
% % % %     range90 = zeros(2, width+1);
% % % %     % 2. Loop A: for each time point, from 0 to 400 (column, left to right)
% % % %     for i = 0:width
% % % %         % 3. Create an array for the value at the point
% % % %         temp_y = y(3:end);
% % % %         values_time = temp_y(find(x(3:end) == i));
% % % %         % 4a. Find the 90%: sort the numbers
% % % %         % 4b. Find the 90%: delete 5% min and 5% max
% % % %         % 4c. Get the range from min and max of existing array
% % % %         values_time_sorted = sort(values_time);
% % % %         values_time_length = length(values_time_sorted);
% % % %         if values_time_length >= 20
% % % %             min_id = values_time_length/20 + 1;
% % % %             max_id = values_time_length*19/20;
% % % %             range90(1,1+i) = values_time_sorted(min_id);
% % % %             range90(2,1+i) = values_time_sorted(max_id);
% % % %         else
% % % %             range90(1,1+i) = 0;
% % % %             range90(2,1+i) = 0;
% % % % %             disp('90% ranges not selected due to simulation size lower than 20.')
% % % %         end
% % % %         % 5. Loop B: for each number in the column (given column, up to down
% % % %         for j = 1:height+1
% % % % %             if((j == 10)&&(i == 86))
% % % % %                 disp('stop here')
% % % % %             end
% % % % %             j/(height+1)
% % % %             temp_dens = [];
% % % %             % 6. Loop C: calculate the inverse of the distance to each data point, and keep adding (the closer, the higher
% % % % %             if((j == 1)||(j == 12))
% % % % %                 disp('stop here')
% % % % %             end
% % % %             for k = 1:values_time_length
% % % %                 fudge = 0;
% % % %                 diff_height = values_time_sorted(k) - (max(y)-min(y))*j/(height+1);
% % % %                 temp_dens = [temp_dens, fudge + 1/(diff_height^(1/2))];
% % % %             end
% % % %             dmap(j,1+i) = sum(abs(temp_dens));
% % % %         end  
% % % %         
% % % %         
% % % % % % % %         % averaging values to eliminate outliers
% % % % % % % %         for j = 1:height+1
% % % % % % % %             if((j == 1)||(j == height+1))
% % % % % % % %             elseif((j == 1+1)||(j == height+1-1))
% % % % % % % % %                 temp_range = [1+i-0 1+i+0];
% % % % % % % %                 temp_range = [j-1 j+1];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+2)||(j == height+1-2))
% % % % % % % % %                 temp_range = [1+i-1 1+i+1];
% % % % % % % %                 temp_range = [j-2 j+2];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+3)||(j == height+1-3))
% % % % % % % % %                 temp_range = [1+i-2 1+i+2];
% % % % % % % %                 temp_range = [j-3 j+3];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+4)||(j == height+1-4))
% % % % % % % % %                 temp_range = [1+i-3 1+i+3];
% % % % % % % %                 temp_range = [j-4 j+4];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+5)||(j == height+1-5))
% % % % % % % % %                 temp_range = [1+i-4 1+i+4];
% % % % % % % %                 temp_range = [j-5 j+5];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+6)||(j == height+1-6))
% % % % % % % % %                 temp_range = [1+i-5 1+i+5];
% % % % % % % %                 temp_range = [j-6 j+6];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+7)||(j == height+1-7))
% % % % % % % % %                 temp_range = [1+i-6 1+i+6];
% % % % % % % %                 temp_range = [j-7 j+7];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+8)||(j == height+1-8))
% % % % % % % % %                 temp_range = [1+i-7 1+i+7];
% % % % % % % %                 temp_range = [j-8 j+8];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j == 1+9)||(j == height+1-9))
% % % % % % % % %                 temp_range = [1+i-8 1+i+8];
% % % % % % % %                 temp_range = [j-9 j+9];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             elseif((j >= 1+9+10)&&(j <= height+1-19))
% % % % % % % %                 temp_range = [j-19 j+19];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             else
% % % % % % % %                 temp_range = [j-10 j+10];
% % % % % % % %                 dmap(j,1+i) = mean(dmap(temp_range,1+i));
% % % % % % % %             end
% % % % % % % %         end  
% % % % % % % %         % 
% % % % % % % % %         if(i == 50)
% % % % % % % % %             disp('stop here')
% % % % % % % % %         end
% % % %     end
    
%     minVal = min(min(dmap));
% %     idxs_temp = find(dmap >= 1000);
%     idxs_temp = find(dmap >= 500);
%     dmap(idxs_temp) = minVal;
    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Hist counts to make the 2D histogram plot
%     X = randn(1000,1);
%     edges = [-5 -4 -2 -1 -0.5 0 0.5 1 2 4 5];
%     N = histcounts(X,edges)
    
    temp_x = x(3:end);
    temp_y = y(3:end);
    temp_edges = [0:1:robust_setup.dpVal_yaxis+1];
    if robust_setup.maxRange_concentration == 0
%         edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_rate;
        % 2022-01-17 edition: to solve the issue of values going above the
        % robust_setup.maxRange_concentration set, and not having to set it
        % manually.
        edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_rate;
        if max(temp_y) >= max(edges)
            robust_setup.maxRange_rate = max(temp_y)*1.1;
            edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_rate;
        end
    else
%         edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_concentration;
        % 2022-01-17 edition: to solve the issue of values going above the
        % robust_setup.maxRange_concentration set, and not having to set it
        % manually.
        edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_concentration;
% %         if(isfield(robust_setup,'specialcase')&&(robust_setup.specialcase ~= 0))
% %             % 
% %             if robust_setup.specialcase == 13 % case of FBP_SS
% %                 robust_setup.maxRange_concentration = 7;
% %                 % 
% %                 indices = find(abs(temp_y)>robust_setup.maxRange_concentration);
% %                 temp_y(indices) = robust_setup.maxRange_concentration;
% %                 % 
% % %                 temp_y = temp_y / max(temp_y) * robust_setup.maxRange_concentration;
% %                 edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_concentration;
% %             end
% %         else
            % 
            if max(temp_y) >= max(edges)
                robust_setup.maxRange_concentration = max(temp_y)*1.1;
                edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_concentration;
            end
% %         end

%         if max(temp_y) >= max(edges)
%             robust_setup.maxRange_concentration = max(temp_y)*1.1;
%             edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_concentration;
%         end

% % % %         % 
% % % %         if(isfield(robust_setup,'forceSSdata')&&(robust_setup.forcedSSdata == 1))
% % % %             edges = temp_edges ./ temp_edges(end) *  robust_setup.maxRange_concentration_forced;
% % % %         end        
        
    end
%     edges = robust_setup.interpSpace;
    dmap = zeros(height+1, width+1);
    range90 = zeros(3, width+1);
    for i = 0:width
        % calculation density map
        temp_idxs = find(temp_x == i);
        temp_dpts = histcounts(temp_y(temp_idxs), edges);
        dmap(:,i+1) = temp_dpts';
        % getting the ranges 0.05-0.95
        range90(1,i+1) = i;
% % % %         range90(2,i+1) = quantile(dmap(:,i+1),0.05);
% % % %         range90(3,i+1) = quantile(dmap(:,i+1),0.95);
        range90(2,i+1) = quantile(temp_y(temp_idxs),0.05);
        range90(3,i+1) = quantile(temp_y(temp_idxs),0.95);
    end
    
%     robust_setup.maxRange_concentration
%     robust_setup.dpVal_yaxis
%     
%     temp1 = [0:1:robust_setup.dpVal_yaxis+1]
%     temp2 = temp1 ./ temp1(end) * 2
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end

