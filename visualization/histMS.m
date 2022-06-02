function [h] = histMS(solutions,setup)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% early size selection
ptot = length(setup.caseStudy.parameters);
for i = 1:10
    j = i^2;
    %qcheck = exist q
    if j >= ptot
        %if(exist(q,'var'))
        exist q var;            % !! pretty ugly this part. Clean up needed. Use 'while' statement...
        k = ans; clear ans
        if k == 1   
        else
            q = sqrt(j);
        end
    end
end

% getting values for the plot
% parray = zeros(length(solutions),ptot); % commented to avoid leaving 0 values
for i = 1:length(solutions)
    if solutions(i).Exitflag == 3
            for j = 1:ptot
                parray(i,j) = solutions(i).X(j);
            end
    end
end

% plotting
h = figure;
for i=1:ptot
    subplot(q,q,i)
    if setup.parEst.MSnonReg == 1
        histogram(parray(:,i),[-3:0.15:3])
        hold on           %to show the limits
        ylims = ylim;
        line([-1 -1],[ylims(1) ylims(2)],'Color','red','LineStyle','--')
        hold on
        line([1 1],[ylims(1) ylims(2)],'Color','red','LineStyle','--')
%         set(gca,'xtick',[],'ytick',[])
    elseif setup.parEst.MSReg == 1
        if(setup.caseStudy.PFK == 1 && (i == 1 || i == 3 || i == 9))
            histogram(parray(:,i),[-3:0.15:3])
            hold on           %to show the limits
            ylims = ylim;
            line([-1 -1],[ylims(1) ylims(2)],'Color','red','LineStyle','--')
            hold on
            line([1 1],[ylims(1) ylims(2)],'Color','red','LineStyle','--')
        else
%             histogram(parray(:,i),[-1:0.1:1])
            histogram(parray(:,i),[-3:0.15:3])
            hold on           %to show the limits
            ylims = ylim;
            line([-1 -1],[ylims(1) ylims(2)],'Color','red','LineStyle','--')
            hold on
            line([1 1],[ylims(1) ylims(2)],'Color','red','LineStyle','--')  
        end
    end
    title(setup.legenda.parameters(setup.caseStudy.parameters(i)))
end

if setup.parEst.MSnonReg == 1
    suptitle(['MultiStart. NO regularized ', setup.legenda.fluxes(setup.caseStudy.fluxes)])
elseif setup.parEst.MSReg == 1
    suptitle(['MultiStart. Regularized ', setup.legenda.fluxes(setup.caseStudy.fluxes)])
end

end

