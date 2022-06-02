function [h] = plotPLA(output,setup)

if setup.parEst.nonreg == 1
    Par     = struct2cell(output.PLAnonReg.pPar);
    Res     = struct2cell(output.PLAnonReg.pRes);
    par     = output.MSparestNonReg.b;
    ptot    = size(struct2cell(output.PLAnonReg.pPar));
elseif setup.parEst.reg == 1
    Par     = struct2cell(output.PLAreg.pPar);
    Res     = struct2cell(output.PLAreg.pRes);
    par     = output.MSparestReg.b;
    ptot    = size(struct2cell(output.PLAreg.pPar));
end

h=figure;
for i = 1:ptot(1)
    plPar = Par{i};
    plRes = Res{i};
    index       = setup.caseStudy.parameters(i);
    parname     = setup.legenda.parameters{index};%name
    
    if ((setup.caseStudy.PYK == 1) || (setup.caseStudy.FBA == 1) || (setup.caseStudy.GPD == 1) || (setup.caseStudy.GLT_HXK == 1))
        subplot(3,3,i)
    elseif ((setup.caseStudy.PFK == 1) || (setup.caseStudy.GAPDH_PGK == 1) || (setup.caseStudy.PDC_ADH == 1))
        subplot(4,4,i)
    else
        subplot(2,2,i)
    end
    plot(plPar, plRes)
    title([parname])
    xlabel('Parameter value')
    ylabel('Fit residual')
    if setup.parEst.nonreg == 1 && setup.parEst.reg == 0
        xlim([-3 3])
    elseif setup.parEst.reg == 1 && setup.parEst.nonreg == 0
        if(setup.caseStudy.PFK == 1 && (i == 1 || i == 3 || i == 9))
            xlim([-3 3])
        else
%             xlim([-1 1])
            xlim([-3 3])
        end
    end
    hold on
    [val,ind] = min(plRes(~ismember(plRes,0)));
    ylim([min(plRes)-0.1 max(plRes)+0.1])
    y1 = ylim;
    line([par(i) par(i)],y1,'Color','magenta','LineStyle',':') % this is the right parameter value
    hold on
    line([plPar(ind) plPar(ind)],y1,'Color','red','LineStyle','--') % this is the current minimum
    hold on
    drawnow()
end

end