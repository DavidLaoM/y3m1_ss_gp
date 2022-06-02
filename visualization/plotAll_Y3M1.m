function [numFigs] = plotAll_Y3M1(simResults,legenda,canelas_SS,data,setup,xAll)
%PLOTALL_Y3M1 This function plots all the variables that give
%information on the sysem simulation. It adds the process variables that
%halp validate a model like this
%   figure 1:   Steady state metabolite profile. All
%   figure 2:   Steady state flux profile. All
%   figure 3:   Dynamic glucose pulse metabolite profile. All
%   figure 4:   Dynamic glucose pulse flux profile. All
%   figure 5:   Dynamic glucose pulse metabolite profile. All. Until 2000s
%   figure 6:   Dynamic glucose pulse flux profile. All. Until 2000s
%   figure 11:  Parameter values. All
%   figure 12:  Parameter values. glycolysis + trehalose + glycerol. Type. Vm and kcat
%   figure 13:  Parameter values. glycolysis + trehalose + glycerol. Type. keq
%   figure 14:  Parameter values. glycolysis + trehalose + glycerol. Type. km
%   figure 15:  Parameter values. glycolysis + trehalose + glycerol. Type. others
%   figure 16:  Parameter values. cofactors and sink reactions
%   figure 21:  Physiology. Steady state profile. All
%   figure 22:  Physiology. Dynamic glucose pulse profile. All
%   figure 23:  Physiology. Steady state profile. Extracellular concentrations
%   figure 24:  Physiology. Dynamic glucose pulse profile. Extracellular concentrations
%   figure 31:  Ratios modules. Steady state profile
%   figure 32:  Ratios modules. Dynamic glucose pulse profile.
%   figure 33:  Imbalance check. Dynamic glucose pulse profile.
%   figure 41:  Metabolite pools. Steady state profile.
%   figure 42:  Metabolite pools. Fluxes related. Steady state profile.
%   figure 43:  Metabolite pools. Dynamic glucose pulse profile.
%   figure 44:  Metabolite pools. Fluxes related. Dynamic glucose pulse profile.
% Structure of the setup options: seup.plotResults = 0
%   setup.plotResults = 0;      No plots are made
%   setup.plotResults = 1;      Simulations. Metabolties and fluxes. Unfocused.
%   setup.plotResults = 2;      Simulations and physiological validation
%   setup.plotResults = 10;     All. Simulations, physiological validation, ratios, parameters, mass balances.
%   setup.plotResults = 20;     Selected cases. When the thing was under development.
%   setup.plotResults = 30;     Old plots. All together.
%   setup.plotResults = 40;     Old style plots. All together.
%   setup.plotResults = 50;     ONLY plots being develped
%   setup.plotResults = 70;     plots made for the BioSB2020 poster
%   setup.plotResults = 80;     temporary tests
%   setup.plotResults = 81;     test for pEfull_r32_difference_11_14
%   setup.plotResults = 99;     standard simulations

% recall data
[npSets,~] = size(xAll);
colorSet = cool(npSets);
recallSimsData_Y3M1;
% default reference parameter array
refParams = exist('setup.refParams');
if refParams == 0 % if 0, first array gets plotted in bold. If-1, none does
    setup.refParams = 1;
end
% background prot color (if exists. If not, default)
if isfield(setup,'plotBackColor')
    backColor = setup.plotBackColor;
else
    backColor = [0.94 0.94 0.94];
end


% plots profiles
if(setup.plotResultsMode == 99)

    % fig1. Steady state metabolite profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Yss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.metabolites.onoff(k) == 1
                errorbar(canelas_SS.ordered.dprof{k},canelas_SS.ordered.metabolites.data{k},canelas_SS.ordered.metabolites.std{k},'r.')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            if k == 25
                xlabel('dilution rate [h^{-1}]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig1. Steady state metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    % fig2. Steady state flux profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state flux profile. All','units','normalized','outerposition',[0 0 1 1])
    
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Vss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.fluxes.onoff(k) == 1
                plot(canelas_SS.ordered.dprof{end},canelas_SS.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            if k == 37
                xlabel('dilution rate [h^{-1}]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end
       
        suptitle('fig2. Steady state flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    %fig3. Dynamic glucose pulse metabolite profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.metabolites.onoff(k) == 1
                plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            xlim([-100 340])
% % % %             xlim([-3000 340])
% % % %             xlim([-100 3400])
            if k == 25
                xlabel('time [s]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig3. Dynamic glucose pulse metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
    
    %fig4. Dynamic glucose pulse flux profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.fluxes.onoff(k) == 1
                plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            xlim([-100 340])
% % % %             xlim([-3000 340])
% % % %             xlim([-100 3400])
            if k == 37
                xlabel('time [s]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end 
        
        % legenda (if necessary)
        if isfield(setup,'namesHits')
%             legendaNames = cell(1,length(npSets));
%             for j = 1:npSets
%                 legendaNames{j} = num2str(j);
%             end
%             lgd = legend(legendaNames);
            lgd = legend(setup.namesHits');
%             lgd.Position = [0.6022 0.0650 0.0305 0.1373];
%             lgd.Position = [1 0.0650 0.0305 0.1373];
            lgd.Position = [0.75 0.0650 0.0305 0.1373];
        end
        
        suptitle('fig4. Dynamic glucose pulse flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    % fig21. Physiology. Steady state profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Physiology. Steady state profile. All','units','normalized','outerposition',[0 0 1 1])
        % % Inputs
        % qS
        subplot(3,6,2)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GLT, 'r+')
        end
        title('q_{S}')
        % mu
        subplot(3,6,7)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), dprofile{j}(1:q(1)), 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), dprofile{j}(1:q(1)), 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.D, 'r+')
        end
        title('growth rate, /mu')
        % % Exchange rates
        % qEtoh
        subplot(3,6,8)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH, 'r+')
        end
        title('q_{Etoh}')
        % qGlycerol
        subplot(3,6,9)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_G3PDH, 'r+')
        end
        title('q_{Glycerol}')
        % qAcet
        subplot(3,6,10)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkACE{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkACE{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.v_PDC - canelas_SS.mtlD.v_ADH), 'r+')
        end
        title('q_{Acet}')
        % qCO2
        subplot(3,6,11)
        for j = 1:npSets
            if j == setup.refParams
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'k','LineWidth',1.2)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qCO2, 'b+')
        end
        title('q_{CO2}')
        % qO2
        subplot(3,6,12)
        for j = 1:npSets
            if j == setup.refParams
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'k','LineWidth',1.2)
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'k','LineWidth',1.2)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'k','LineWidth',1.2)
            end
        hold on
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qO2, 'b+')
        end
        title('q_{O2}')
% % % %      
        % calc qmito and qATPase
        for k = 27
            qATPase = zeros(1,length(dprofile{1}));
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{1}.ss2can{o}(end,k);
                qmito(o) = Vss{1}.ss2can{o}(end,k+1);
            end
        end
        % qmito
        subplot(3,6,5)
        for j = 1:npSets
            % calc qmito and qATPase
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qmito(o) = Vss{j}.ss2can{o}(end,k+1);
            end            
            % sims            
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), qmito, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), qmito, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D,(exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]'/1, 'b+')
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.sinkPYR * 3 + 1 * canelas_SS.mtlD.v_PDC)*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]/1, 'b+')
        end
        title('q_{mito}')
        
        % qATPase
        subplot(3,6,6)
        for j = 1:npSets
% % % %                
            qATPase = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{j}.ss2can{o}(end,k);
            end            
% % % %             
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), qATPase, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), qATPase, 'color', colorSet(j,:))
            end
            hold on
%             plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/1.7/3.6, 'b+')
            plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'b+')
            hold on
            plot(canelas_SS.mtlD.D,(40*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'g+')
        end
        title('q_{ATPase}')
% % % %         
        % Yxs
        subplot(3,6,13)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, dprofile{j}./exp_v_GLT{j}', 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, dprofile{j}./exp_v_GLT{j}', 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.D./canelas_SS.mtlD.v_GLT), 'r+')
        end        
        title('Y_{XS}')
        % Yps
        subplot(3,6,14)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.v_ADH./canelas_SS.mtlD.v_GLT), 'r+')
        end    
        title('Y_{PS}')
        % % Ratios
        % POratio: recalculate back
        subplot(3,6,15)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_mito{j}./(exp_v_mito{j}/1.1), 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_mito{j}./(exp_v_mito{j}/1.1), 'color', colorSet(j,:))
            end
        hold on
        end    
        title('PO_{ratio}')
        % RQratio
        subplot(3,6,16)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'color', colorSet(j,:))
            end
        hold on
        end    
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qCO2./canelas_SS.mtlD.qO2, 'b+');
        title('RQ_{ratio}')
        % redox state
        subplot(3,6,17)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_NAD{j}./exp_NADH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_NAD{j}./exp_NADH{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.NAD_NADHratio, 'r+')
        end    
        title('NAD:NADH_{ratio}')
        
        subplot(3,6,18)
        for j = 1:npSets
            sumAXP = exp_ATP{j} + exp_ADP{j} + exp_AMP{j};
            if j == setup.refParams
                plot(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'k','LineWidth',1.2)
            else
                plot(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP/2) ./ (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP + canelas_SS.mtlD.AMP), 'r+')
        end
        title('Energy charge: (ATP + ADP/2) / sumAXP')
        
        % set(gcf,'color','w'); suptitle('fig21. Physiology. Steady state profile. All')
        set(gcf,'color',backColor);
    end
    end
        
    
end

% plots profiles
if((setup.plotResultsMode == 1)||(setup.plotResultsMode == 2)||(setup.plotResultsMode == 10)||(setup.plotResultsMode == 30))

    % fig1. Steady state metabolite profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Yss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.metabolites.onoff(k) == 1
                errorbar(canelas_SS.ordered.dprof{k},canelas_SS.ordered.metabolites.data{k},canelas_SS.ordered.metabolites.std{k},'r.')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            if k == 25
                xlabel('dilution rate [h^{-1}]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig1. Steady state metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    % fig2. Steady state flux profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state flux profile. All','units','normalized','outerposition',[0 0 1 1])
    
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Vss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.fluxes.onoff(k) == 1
                plot(canelas_SS.ordered.dprof{end},canelas_SS.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            if k == 37
                xlabel('dilution rate [h^{-1}]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end
       
        suptitle('fig2. Steady state flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    %fig3. Dynamic glucose pulse metabolite profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.metabolites.onoff(k) == 1
                plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            xlim([-100 340])
% % % %             xlim([-3000 340])
% % % %             xlim([-100 3400])
            if k == 25
                xlabel('time [s]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig3. Dynamic glucose pulse metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
    
    %fig4. Dynamic glucose pulse flux profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.fluxes.onoff(k) == 1
                plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            xlim([-100 340])
% % % %             xlim([-3000 340])
% % % %             xlim([-100 3400])
            if k == 37
                xlabel('time [s]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end 
        
        % legenda (if necessary)
        if isfield(setup,'namesHits')
%             legendaNames = cell(1,length(npSets));
%             for j = 1:npSets
%                 legendaNames{j} = num2str(j);
%             end
%             lgd = legend(legendaNames);
            lgd = legend(setup.namesHits');
%             lgd.Position = [0.6022 0.0650 0.0305 0.1373];
%             lgd.Position = [1 0.0650 0.0305 0.1373];
            lgd.Position = [0.75 0.0650 0.0305 0.1373];
        end
        
        suptitle('fig4. Dynamic glucose pulse flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
    
    %fig5. Dynamic glucose pulse metabolite profile. All. Until 2000s
    for i = 2
    if setup.runGSvanHeerden == 2
        figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.metabolites.onoff(k) == 1
                plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            xlim([-100 2000])
            if k == 25
                xlabel('time [s]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig5. Dynamic glucose pulse metabolite profile. All. Long term')
        % set(gcf,'color','w');
    end
    end
    
    %fig6. Dynamic glucose pulse flux profile. All. Until 2000s
    for i = 2
    if setup.runGSvanHeerden == 2
        figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.fluxes.onoff(k) == 1
                plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            xlim([-100 2000])
            if k == 37
                xlabel('time [s]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end       
        suptitle('fig6. Dynamic glucose pulse flux profile. All. Long term.')
        % set(gcf,'color','w');
    end
    end

end


% plots parameters
if(setup.plotResultsMode == 10)
    x = xAll(setup.refParams,:);
    
    % fig11. Parameter values. All
    for i = 1

        figure('name','Parameter values. All','units','normalized','outerposition',[0 0 1 1])
        subplot(1,2,1)
        for j = 1:length(x)
            load('customColormap.mat','customColormap')
            color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,10.^x(j),'o','MarkerFaceColor',color,'MarkerEdgeColor',color)
            hold on
        end
        xlabel('Parameters per ID')
        ylabel('Ration P_i_n_v_i_v_o:P_i_n_v_i_t_r_o')
        ylim([1E-3 1000])

        [p] = convertX2P(x);
        subplot(1,2,2)
        for j = 1:length(p)
            load('customColormap.mat','customColormap')
    %         color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,p(j),'o','color','black')
            hold on
        end
        xlabel('Parameters per ID')
        ylabel('Parameter values')
    %     ylim([1E-3 1000])

        suptitle('fig11. Parameter values. All')
        % set(gcf,'color','w');
    end

    % fig12. Parameter values. glycolysis+trehalose+glycerol. Type. Vm and kcat
    for i = 1
        figure('name','Parameter values. glycolysis+trehalose+glycerol. Type. Vm and kcat','units','normalized','outerposition',[0 0 1 1])
        % vms and kcats
            % 37 34 56 60 15 82 70 19 77 27 66 104 42 10 107 86 86 126 120 123
        leg_vmkc = [37 34 56 60 15 82 70 19 77 27 66 104 42 10 107 86 126 120 123];
        subplot(1,2,1)
        for j = 1:length(leg_vmkc)
            color = customColormap(int8(1+127*(max(abs(x(leg_vmkc(j))))/max(abs(x)))),:);
            semilogy(j,10.^x(leg_vmkc(j)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
            hold on
        end
        ylabel('Ration P_i_n_v_i_v_o:P_i_n_v_i_t_r_o')
        ylim([1E-3 1000])
        htemp = gca;
        xlim([0 htemp.XLim(2)+1])
        line([0 htemp.XLim(2)],[1E-1 1E-1],'Color','black','LineStyle','--')
        line([0 htemp.XLim(2)],[1E1 1E1],'Color','black','LineStyle','--')
        for k = 1:length(leg_vmkc)
%             ptemp_val = p(leg_vmkc(k));
            ptemp_xval = x(leg_vmkc(k));
            ptemp_idx = k;
            ptemp_name = legenda.eName{leg_vmkc(k)};
            text(0.2+ptemp_idx,10^ptemp_xval,ptemp_name)
        end
        
        subplot(1,2,2)
        for j = 1:length(leg_vmkc)
            load('customColormap.mat','customColormap')
    %         color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,p(leg_vmkc(j)),'o','color','black','MarkerFaceColor','black')
            hold on
        end
    %     xlabel('Parameters per ID')
        ylabel('Parameter values')
    %     ylim([1E-3 1000])   
        htemp = gca;
        xlim([0 htemp.XLim(2)+1])
        for k = 1:length(leg_vmkc)
%             ptemp_val = p(leg_vmkc(k));
            ptemp_pval = p(leg_vmkc(k));
            ptemp_idx = k;
            ptemp_name = legenda.eName{leg_vmkc(k)};
            text(0.2+ptemp_idx,ptemp_pval,ptemp_name)
        end

        set(gcf,'color','w'); suptitle('fig12. Parameter values. glycolysis+trehalose+glycerol. Type. Vm and kcat')

    end

    % fig13. Parameter values. glycolysis+trehalose+glycerol. Type. keq
    for i = 1
        figure('name','Parameter values. glycolysis+trehalose+glycerol. Type. keq','units','normalized','outerposition',[0 0 1 1])
        % keqs
            % 30 57 12 80 69 17 21 61 99 1 83 83
        leg_keq = [30 57 12 80 69 17 21 61 99 1 83];
        subplot(1,2,1)
        for j = 1:length(leg_keq)
            color = customColormap(int8(1+127*(max(abs(x(leg_keq(j))))/max(abs(x)))),:);
            semilogy(j,10.^x(leg_keq(j)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color)
            hold on
        end
        ylabel('Ration P_i_n_v_i_v_o:P_i_n_v_i_t_r_o')
        ylim([1E-3 1000])
        htemp = gca;
        xlim([0 htemp.XLim(2)+1])
        line([0 htemp.XLim(2)],[1E-1 1E-1],'Color','black','LineStyle','--')
        line([0 htemp.XLim(2)],[1E1 1E1],'Color','black','LineStyle','--')
        for k = 1:length(leg_keq)
%             ptemp_val = p(leg_keq(k));
            ptemp_xval = x(leg_keq(k));
            ptemp_idx = k;
            ptemp_name = legenda.eName{leg_keq(k)};
            text(0.2+ptemp_idx,10^ptemp_xval,ptemp_name)
        end


        subplot(1,2,2)
        for j = 1:length(leg_keq)
            load('customColormap.mat','customColormap')
    %         color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,p(leg_keq(j)),'o','color','black','MarkerFaceColor','black')
            hold on
        end
    %     xlabel('Parameters per ID')
        ylabel('Parameter values')
    %     ylim([1E-3 1000])  
        htemp = gca;
        xlim([0 htemp.XLim(2)+1])
        for k = 1:length(leg_keq)
%             ptemp_val = p(leg_keq(k));
            ptemp_pval = p(leg_keq(k));
            ptemp_idx = k;
            ptemp_name = legenda.eName{leg_keq(k)};
            text(0.2+ptemp_idx,ptemp_pval,ptemp_name)
        end 


        set(gcf,'color','w'); suptitle('fig13. Parameter values. glycolysis+trehalose+glycerol. Type. keq')
    end

    % fig14. Parameter values. glycolysis+trehalose+glycerol. Type. km
    for i = 1
        figure('name','Parameter values. glycolysis+trehalose+glycerol. Type. km','units','normalized','outerposition',[0 0 1 1])
        % kms
            % 35 36 28 29 31 32 33 48 49 50 51 52 58 59 11 13 14 79 81 67 68 16 18
            % 71 72 73 74 22 23 24 25 26 62 63 64 65 96 97 98 100 101 102 103 39 40
            % 6 7 8 9 105 106 84 85 108 95 84 85 124 125 127 128 119 121 122
        leg_km = [...
                    35 36 28 29 31 32 33 48 49 50 51 52 58 59 11 13 14 79 81 67 68 16 18,...
                    71 72 73 74 22 23 24 25 26 62 63 64 65 96 97 98 100 101 102 103 39 40,...
                    6 7 8 9 105 106 84 85 108 84 85 124 125 127 128 119 121 122,...
                    ];
        subplot(1,2,1)
        for j = 1:length(leg_km)
            color = customColormap(int8(1+127*(max(abs(x(leg_km(j))))/max(abs(x)))),:);
            semilogy(j,10.^x(leg_km(j)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color)
            hold on
        end
        ylabel('Ration P_i_n_v_i_v_o:P_i_n_v_i_t_r_o')
        ylim([1E-3 1000])
        htemp = gca;
        xlim([0 htemp.XLim(2)+1])
        line([0 htemp.XLim(2)],[1E-1 1E-1],'Color','black','LineStyle','--')
        line([0 htemp.XLim(2)],[1E1 1E1],'Color','black','LineStyle','--')
        for k = 1:length(leg_km)
%             ptemp_val = p(leg_km(k));
            ptemp_xval = x(leg_km(k));
            ptemp_idx = k;
            ptemp_name = legenda.eName{leg_km(k)};
            text(0.2+ptemp_idx,10^ptemp_xval,ptemp_name)
        end

        subplot(1,2,2)
        for j = 1:length(leg_km)
            load('customColormap.mat','customColormap')
    %         color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,p(leg_km(j)),'o','color','black','MarkerFaceColor','black')
            hold on
        end
    %     xlabel('Parameters per ID')
        ylabel('Parameter values')
        htemp = gca;
        xlim([0 htemp.XLim(2)+1])
        for k = 1:length(leg_km)
%             ptemp_val = p(leg_km(k));
            ptemp_pval = p(leg_km(k));
            ptemp_idx = k;
            ptemp_name = legenda.eName{leg_km(k)};
            text(0.2+ptemp_idx,ptemp_pval,ptemp_name)
        end

        set(gcf,'color','w'); suptitle('fig14. Parameter values. glycolysis+trehalose+glycerol. Type. km')
    end

    % fig15. Parameter values. glycolysis+trehalose+glycerol. Type. others
    for i = 1
        figure('name','Parameter values. glycolysis+trehalose+glycerol. Type. others','units','normalized','outerposition',[0 0 1 1])
        % just check others
            % 43 44 45 46 47 53 54 55 87 75 76 41 2 3 4 5
        leg_others = [43 44 45 46 47 53 54 55 87 75 76 41 2 3 4 5];

        subplot(1,2,1)
        for j = 1:length(leg_others)
            color = customColormap(int8(1+127*(max(abs(x(leg_others(j))))/max(abs(x)))),:);
            semilogy(j,10.^x(leg_others(j)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color)
            hold on
        end
        ylabel('Ration P_i_n_v_i_v_o:P_i_n_v_i_t_r_o')
        ylim([1E-3 1000])

        subplot(1,2,2)
        for j = 1:length(leg_others)
            load('customColormap.mat','customColormap')
    %         color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,p(leg_others(j)),'o','color','black')
            hold on
        end
    %     xlabel('Parameters per ID')
        ylabel('Parameter values') 

        % set(gcf,'color','w'); suptitle('fig15. Parameter values. glycolysis+trehalose+glycerol. Type. others')
    end

    % fig16. Parameter values. cofactors
    for i = 1
        figure('name','Parameter values. cofactors','units','normalized','outerposition',[0 0 1 1])
        % just wrte all from the cofactors
            % 88 94 89 90 91 92 93 109 110 111 129 112 113 114 115 116 117 118 37
            % 38 130 131
        leg_cof = [
                    88 94 89 90 91 92 93 109 110 111 129 112 113 114 115 116 117 118 37,...
                    38 130 131,...
                    ];
        subplot(1,2,1)
        for j = 1:length(leg_cof)
            if length(x) == 128
                x(129:131) = zeros;
            end
                color = customColormap(int8(1+127*(max(abs(x(leg_cof(j))))/max(abs(x)))),:);
                semilogy(j,10.^x(leg_cof(j)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color)
                hold on
        end
        ylabel('Ration P_i_n_v_i_v_o:P_i_n_v_i_t_r_o')
        ylim([1E-3 1000])

        subplot(1,2,2)
        for j = 1:length(leg_cof)
            load('customColormap.mat','customColormap')
    %         color = customColormap(int8(1+127*(max(abs(x(j)))/max(abs(x)))),:);
            semilogy(j,p(leg_cof(j)),'o','color','black')
            hold on
        end
    %     xlabel('Parameters per ID')
        ylabel('Parameter values') 


        % set(gcf,'color','w'); suptitle('fig16. Parameter values. cofactors ')
    end

end


% plots physiology
if((setup.plotResultsMode == 2)||(setup.plotResultsMode == 10)||(setup.plotResultsMode == 20)||(setup.plotResultsMode == 30))

    % fig21. Physiology. Steady state profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Physiology. Steady state profile. All','units','normalized','outerposition',[0 0 1 1])
        % % Inputs
        % qS
        subplot(3,6,2)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GLT, 'r+')
        end
        title('q_{S}')
        % mu
        subplot(3,6,7)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), dprofile{j}(1:q(1)), 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), dprofile{j}(1:q(1)), 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.D, 'r+')
        end
        title('growth rate, /mu')
        % % Exchange rates
        % qEtoh
        subplot(3,6,8)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH, 'r+')
        end
        title('q_{Etoh}')
        % qGlycerol
        subplot(3,6,9)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_G3PDH, 'r+')
        end
        title('q_{Glycerol}')
        % qAcet
        subplot(3,6,10)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkACE{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkACE{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.v_PDC - canelas_SS.mtlD.v_ADH), 'r+')
        end
        title('q_{Acet}')
        % qCO2
        subplot(3,6,11)
        for j = 1:npSets
            if j == setup.refParams
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'k','LineWidth',1.2)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qCO2, 'b+')
        end
        title('q_{CO2}')
        % qO2
        subplot(3,6,12)
        for j = 1:npSets
            if j == setup.refParams
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'k','LineWidth',1.2)
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'k','LineWidth',1.2)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'k','LineWidth',1.2)
            end
        hold on
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qO2, 'b+')
        end
        title('q_{O2}')
% % % %      
        % calc qmito and qATPase
        for k = 27
            qATPase = zeros(1,length(dprofile{1}));
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{1}.ss2can{o}(end,k);
                qmito(o) = Vss{1}.ss2can{o}(end,k+1);
            end
        end
        % qmito
        subplot(3,6,5)
        for j = 1:npSets
            % calc qmito and qATPase
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qmito(o) = Vss{j}.ss2can{o}(end,k+1);
            end            
            % sims            
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), qmito, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), qmito, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D,(exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]'/1, 'b+')
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.sinkPYR * 3 + 1 * canelas_SS.mtlD.v_PDC)*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]/1, 'b+')
        end
        title('q_{mito}')
        
        % qATPase
        subplot(3,6,6)
        for j = 1:npSets
% % % %                
            qATPase = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{j}.ss2can{o}(end,k);
            end            
% % % %             
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), qATPase, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), qATPase, 'color', colorSet(j,:))
            end
            hold on
%             plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/1.7/3.6, 'b+')
            plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'b+')
            hold on
            plot(canelas_SS.mtlD.D,(40*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'g+')
        end
        title('q_{ATPase}')
% % % %         
        % Yxs
        subplot(3,6,13)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, dprofile{j}./exp_v_GLT{j}', 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, dprofile{j}./exp_v_GLT{j}', 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.D./canelas_SS.mtlD.v_GLT), 'r+')
        end        
        title('Y_{XS}')
        % Yps
        subplot(3,6,14)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.v_ADH./canelas_SS.mtlD.v_GLT), 'r+')
        end    
        title('Y_{PS}')
        % % Ratios
        % POratio: recalculate back
        subplot(3,6,15)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_mito{j}./(exp_v_mito{j}/1.1), 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_mito{j}./(exp_v_mito{j}/1.1), 'color', colorSet(j,:))
            end
        hold on
        end    
        title('PO_{ratio}')
        % RQratio
        subplot(3,6,16)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'color', colorSet(j,:))
            end
        hold on
        end    
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qCO2./canelas_SS.mtlD.qO2, 'b+');
        title('RQ_{ratio}')
        % redox state
        subplot(3,6,17)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_NAD{j}./exp_NADH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_NAD{j}./exp_NADH{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.NAD_NADHratio, 'r+')
        end    
        title('NAD:NADH_{ratio}')
        
        subplot(3,6,18)
        for j = 1:npSets
            sumAXP = exp_ATP{j} + exp_ADP{j} + exp_AMP{j};
            if j == setup.refParams
                plot(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'k','LineWidth',1.2)
            else
                plot(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP/2) ./ (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP + canelas_SS.mtlD.AMP), 'r+')
        end
        title('Energy charge: (ATP + ADP/2) / sumAXP')
        
        % set(gcf,'color','w'); suptitle('fig21. Physiology. Steady state profile. All')
        set(gcf,'color',backColor);
    end
    end
    
    % fig22. Physiology. Dynamic glucose pulse profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Physiology. Dynamic glucose pulse profile. All','units','normalized','outerposition',[0 0 1 1])
        % % Inputs
        % qS
        subplot(3,6,2)
%         plot(Tgs,Vgs(:,1),'-','color','black'), xlim([-100 340])
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,1), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,1), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('q_{S}')
        % mu
        subplot(3,6,7)
%         plot(Tgs,dprofile(3)*ones(size(Tgs)),'-','color','black'), xlim([-100 340])
        if setup.runSScanelas == 1
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, dprofile{1}(3)*ones(size(Tgs{j})), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, dprofile{1}(3)*ones(size(Tgs{j})), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        end
        title('growth rate, /mu')
        % % Exchange rates
        % qEtoh
        subplot(3,6,8)
%         plot(Tgs,Vgs(:,23),'-','color','black'), xlim([-100 340])
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,23), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,23), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('q_{Etoh}')
        % qGlycerol
        subplot(3,6,9)
%         plot(Tgs,Vgs(:,24),'-','color','black'), xlim([-100 340])
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,24), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,24), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('q_{Glycerol}')
        % qAcet
        subplot(3,6,10)
%         plot(Tgs,-Vgs(:,40),'-','color','black'), xlim([-100 340])
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, -Vgs{j}(:,40), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, -Vgs{j}(:,40), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('q_{Acet}')
        % qCO2
        subplot(3,6,11)
%         plot(Tgs,-Vgs(:,39) * 2 + 1 * Vgs(:,13),'-','color','black'), xlim([-100 340]) % vsinkPYR*2+vPDC*1
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,39) * 3 + 1 * Vgs{j}(:,13), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,39) * 3 + 1 * Vgs{j}(:,13), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('q_{CO2}')
        % qO2
        subplot(3,6,12)
%         plot(Tgs,Vgs(:,28)/1.1,'-','color','black'), xlim([-100 340]) % vmito/1.1
        for j = 1:npSets
            if j == setup.refParams
%                 plot(Tgs{j}, Vgs{j}(:,28)/1.1, 'k','LineWidth',1.2), xlim([-100 340])
                plot(Tgs{j}, (Vgs{j}(:,39) * 3 + 1 * Vgs{j}(:,13))/1.1, 'k','LineWidth',1.2), xlim([-100 340])
            else
%                 plot(Tgs{j}, Vgs{j}(:,28)/1.1, 'color', colorSet(j,:)), xlim([-100 340])
                plot(Tgs{j}, (Vgs{j}(:,39) * 3 + 1 * Vgs{j}(:,13))/1.1, 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('q_{O2}')
        % Yxs
        subplot(3,6,13)
%         plot(Tgs,dprofile(3)./Vgs(:,1)','-','color','black'), xlim([-100 340]) % mu/qS
        if setup.runSScanelas == 1
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, dprofile{1}(3)./Vgs{j}(:,1), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, dprofile{1}(3)./Vgs{j}(:,1), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        end
        title('Y_{XS}')
        % Yps
        subplot(3,6,14)
%         plot(Tgs,Vgs(:,23)./Vgs(:,1),'-','color','black'), xlim([-100 340]) % qEtoh/qS
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,23)./Vgs{j}(:,1), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,23)./Vgs{j}(:,1), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('Y_{PS}')
        % % Ratios
        % POratio: recalculate back
        subplot(3,6,15)
%         plot(Tgs,Vgs(:,28)./(Vgs(:,28)/1.1),'-','color','black'), xlim([-100 340]) % vmito/qO2 (it will just be a check as of now)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,28)./(Vgs{j}(:,28)/1.1), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,28)./(Vgs{j}(:,28)/1.1), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('PO_{ratio}')
        % RQratio
        subplot(3,6,16)
%         plot(Tgs,(-Vgs(:,39) * 2 + 1 * Vgs(:,13)) ./ (Vgs(:,28)/1.1),'-','color','black'), xlim([-100 340]) % qCO2/qO2
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, (-Vgs{j}(:,39) * 2 + 1 * Vgs{j}(:,13)) ./ (Vgs{j}(:,28)/1.1), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, (-Vgs{j}(:,39) * 2 + 1 * Vgs{j}(:,13)) ./ (Vgs{j}(:,28)/1.1), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('RQ_{ratio}')
        % redox state
        subplot(3,6,17)
%         plot(Tgs,Ygs(:,7)./Ygs(:,8),'-','color','black'), xlim([-100 340]) % NAD/NADH
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,7)./Ygs{j}(:,8), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Ygs{j}(:,7)./Ygs{j}(:,8), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        title('NAD:NADH_{ratio}')
        
        subplot(3,6,18)
        for j = 1:npSets
            sumAXP = Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16);
            if j == setup.refParams
                plot(Tgs{j},(Ygs{j}(:,9) + Ygs{j}(:,15)/2)./sumAXP,'k','LineWidth',1.2) 
            else
                plot(Tgs{j},(Ygs{j}(:,9) + Ygs{j}(:,15)/2)./sumAXP,'color', colorSet(j,:)), 
            end
            hold on
            sumAXP = data.nucleotides(:,2) + data.nucleotides(:,3) + data.nucleotides(:,4); 
            plot(data.time_nucleotides,(data.nucleotides(:,4) + data.nucleotides(:,3)/2)./sumAXP,'r+')
        end
        xlim([-100 340]), 
        title('Energy charge: (ATP + ADP/2) / sumAXP')
        
        % set(gcf,'color','w'); suptitle('fig22. Physiology. Dynamic glucose pulse profile. All')
        set(gcf,'color',backColor);
    end
    end

    % fig23. Physiology. Steady state profile. Extracellular concentrations
    % figure('name','Physiology. Steady state profile. Extracellular concentrations','units','normalized','outerposition',[0 0 1 1])
    % set(gcf,'color','w'); suptitle('fig23. Physiology. Steady state profile. Extracellular concentrations')


    % fig24. Physiology. Dynamic glucose pulse profile. Extracellular concentrations
    % figure('name','Physiology. Dynamic glucose pulse profile. Extracellular concentrations','units','normalized','outerposition',[0 0 1 1])
    % set(gcf,'color','w'); suptitle('fig24. Physiology. Dynamic glucose pulse profile. Extracellular concentrations')

end


% plots ratios, imbalance and metabolite pools
if((setup.plotResultsMode == 10)||(setup.plotResultsMode == 20))

    % fig31. Ratios modules. Steady state profile
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Ratios modules. Steady state profile','units','normalized','outerposition',[0 0 1 1])

        subplot(3,4,2) %UG/LG, FBA/GAPDH
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_FBA{j}./exp_v_GAPDH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_FBA{j}./exp_v_GAPDH{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_FBA./canelas_SS.mtlD.v_GAPDH, 'r+')
        end        
        title('UG/LG, FBA/GAPDH')

        subplot(3,4,3) %PDC/ADH
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_PDC{j}./exp_v_ADH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_PDC{j}./exp_v_ADH{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PDC./canelas_SS.mtlD.v_ADH, 'r+')
        end
        title('PDC/ADH')

        subplot(3,4,5) %sinkG6P/vPGI
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkG6P{j}./exp_v_PGI{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkG6P{j}./exp_v_PGI{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.sinkG6P./canelas_SS.mtlD.v_PGI, 'r+')
        end
        title('sinkG6P/vPGI')

        subplot(3,4,6) %sinkF6P/vPFK
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkF6P{j}./exp_v_PFK{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkF6P{j}./exp_v_PFK{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, -canelas_SS.mtlD.sinkF6P./canelas_SS.mtlD.v_PFK, 'r+')
        end
        title('sinkF6P/vPFK')

        subplot(3,4,7) %sinkGAP/vGAPDH
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkGAP{j}./exp_v_GAPDH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkGAP{j}./exp_v_GAPDH{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, -canelas_SS.mtlD.sinkGAP./canelas_SS.mtlD.v_GAPDH, 'r+')
        end        
        title('sinkGAP/vGAPDH')

        subplot(3,4,8) %sinkP3G/vPGM
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkP3G{j}./exp_v_PGM{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkP3G{j}./exp_v_PGM{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.sinkP3G./canelas_SS.mtlD.v_PGM, 'r+')
        end  
        title('sinkP3G/vPGM')

        subplot(3,4,9) %sinkPEP/vPYK
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkPEP{j}./exp_v_PYK{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkPEP{j}./exp_v_PYK{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.sinkPEP./canelas_SS.mtlD.v_PYK, 'r+')
        end          
        title('sinkPEP/vPYK')    

        subplot(3,4,10) %sinkPYR/vPDC
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkPYR{j}./exp_v_PDC{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkPYR{j}./exp_v_PDC{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.sinkPYR./canelas_SS.mtlD.v_PDC, 'r+')
        end
        title('sinkPYR/vPDC')

        subplot(3,4,11) %sinkACE/vADH
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}, exp_v_sinkACE{j}./exp_v_ADH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_sinkACE{j}./exp_v_ADH{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.sinkACE./canelas_SS.mtlD.v_ADH, 'r+')
        end
        title('sink/after ACE')

        % set(gcf,'color','w'); suptitle('fig31. Ratios modules. Steady state profile')
    end
    end
    
    % fig32. Ratios modules. Dynamic glucose pulse profile
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Ratios modules. Dynamic glucose pulse profile','units','normalized','outerposition',[0 0 1 1])

        subplot(3,4,2) %UG/LG, FBA/GAPDH
%         plot(Tgs,Vgs(:,5)./Vgs(:,8),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,5)./Vgs{j}(:,8), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,5)./Vgs{j}(:,8), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('UG/LG, FBA/GAPDH')

        subplot(3,4,3) %PDC/ADH
%         plot(Tgs,Vgs(:,13)./Vgs(:,41),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,13)./Vgs{j}(:,41), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,13)./Vgs{j}(:,41), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('PDC/ADH')

        subplot(3,4,5) %sinkG6P/vPGI
%         plot(Tgs,Vgs(:,34)./Vgs(:,3),'-','color','black')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,34)./Vgs{j}(:,3), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,34)./Vgs{j}(:,3), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end 
        xlim([-100 340])
        title('sinkG6P/vPGI')

        subplot(3,4,6) %sinkF6P/vPFK
%         plot(Tgs,Vgs(:,35)./Vgs(:,4),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,35)./Vgs{j}(:,4), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,35)./Vgs{j}(:,4), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('sinkF6P/vPFK')

        subplot(3,4,7) %sinkGAP/vGAPDH
%         plot(Tgs,Vgs(:,36)./Vgs(:,8),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,36)./Vgs{j}(:,8), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,36)./Vgs{j}(:,8), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('sinkGAP/vGAPDH')

        subplot(3,4,8) %sinkP3G/vPGM
%         plot(Tgs,Vgs(:,37)./Vgs(:,10),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,37)./Vgs{j}(:,10), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,37)./Vgs{j}(:,10), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('sinkP3G/vPGM')

        subplot(3,4,9) %sinkPEP/vPYK
%         plot(Tgs,Vgs(:,38)./Vgs(:,12),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,38)./Vgs{j}(:,12), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,38)./Vgs{j}(:,12), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('sinkPYR/vPDC')    

        subplot(3,4,10) %sinkPYR/vPDC
%         plot(Tgs,Vgs(:,39)./Vgs(:,13),'-','color','black')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,39)./Vgs{j}(:,13), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,39)./Vgs{j}(:,13), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end 
        xlim([-100 340])
        title('sinkPYR/vPDC')

        subplot(3,4,11) %sinkACE/vADH
%         plot(Tgs,Vgs(:,40)./Vgs(:,41),'-','color','black') 
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,40)./Vgs{j}(:,41), 'k','LineWidth',1.2), xlim([-100 340])
            else
                plot(Tgs{j}, Vgs{j}(:,40)./Vgs{j}(:,41), 'color', colorSet(j,:)), xlim([-100 340])
            end
            hold on
        end
        xlim([-100 340])
        title('sink/after ACE')

        % set(gcf,'color','w'); suptitle('fig32. Ratios modules. Dynamic glucose pulse profile')
    end
    end

    % fig33. Imbalance check. Dynamic glucose pulse profile
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Imbalance check. Dynamic glucose pulse profile','units','normalized','outerposition',[0 0 1 1])

        % metabolites
            %G6P
            subplot(3,7,1)
%             plot(Tgs,Ygs(:,5),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,5), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Ygs{j}(:,5), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('G6P')
            %F6P
            subplot(3,7,2)
%             plot(Tgs,Ygs(:,4),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,4), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Ygs{j}(:,4), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('F6P')
            %FBP
            subplot(3,7,3)
%             plot(Tgs,Ygs(:,3),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,3), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Ygs{j}(:,3), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('FBP')
            %GAP
            subplot(3,7,4)
%             plot(Tgs,Ygs(:,14),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,14), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Ygs{j}(:,14), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('GAP')
            %DHAP
            subplot(3,7,5)
%             plot(Tgs,Ygs(:,17),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,17), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Ygs{j}(:,17), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('DHAP')
            %BPG
            subplot(3,7,6)
%             plot(Tgs,Ygs(:,2),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,2), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Ygs{j}(:,2), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('BPG')

        % fluxes
            %vHK
            subplot(3,7,8)
%             plot(Tgs,Vgs(:,2),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,2), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,2), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vHK')
            %vPGI
            subplot(3,7,9)
%             plot(Tgs,Vgs(:,3),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,3), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,3), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vPGI')
            %vPFK
            subplot(3,7,10)
%             plot(Tgs,Vgs(:,4),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,4), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,4), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vPFK')
            %vFBA
            subplot(3,7,11)
%             plot(Tgs,Vgs(:,5),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,5), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,5), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vFBA')
            %vTPI
            subplot(3,7,12)
%             plot(Tgs,Vgs(:,7),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,7), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,7), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vTPI')
            %vG3PDH
            subplot(3,7,13)
%             plot(Tgs,Vgs(:,6),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,6), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,6), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vG3PDH')
            %vGAPDH
            subplot(3,7,14)
%             plot(Tgs,Vgs(:,8),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,8), 'k','LineWidth',1.2), 
                else
                    plot(Tgs{j}, Vgs{j}(:,8), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]),title('vGAPDH')

        % set(gcf,'color','w'); suptitle('fig33. Imbalance check. Dynamic glucose pulse profile')
    end
    end

    % fig41. Metabolite pools. Steady state profile
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Metabolite pools. Steady state profile','units','normalized','outerposition',[0 0 1 1])
        % C-structures
            % sum C-structures
            % glc,g6p,f6p,fbp,gap,dhap,g3p,glyc,g1p,t6p,tre,bpg,3pg,2pg,pep,pyr,ace,etoh
            subplot(5,5,2)
            for j = 1:npSets
                if j == setup.refParams
                    sumCstruct = exp_GLCi{j} + exp_G6P{j} + exp_F6P{j} + exp_FBP{j} + exp_GAP{j} + exp_DHAP{j} + exp_G3P{j} + exp_G1P{j} + exp_T6P{j} + exp_BPG{j} + exp_P3G{j} + exp_P2G{j} + exp_PEP{j} + exp_PYR{j} + exp_ETOH{j} + exp_ACE{j} + exp_GLYC{j} + exp_TRE{j};
                    plot(dprofile{j}, sumCstruct,'-o','color','black') 
                else
                    sumCstruct = exp_GLCi{j} + exp_G6P{j} + exp_F6P{j} + exp_FBP{j} + exp_GAP{j} + exp_DHAP{j} + exp_G3P{j} + exp_G1P{j} + exp_T6P{j} + exp_BPG{j} + exp_P3G{j} + exp_P2G{j} + exp_PEP{j} + exp_PYR{j} + exp_ETOH{j} + exp_ACE{j} + exp_GLYC{j} + exp_TRE{j};
                    plot(dprofile{j}, sumCstruct, 'color', colorSet(j,:))
                end
            hold on
            end
            title('sum C-structures')
        % Trehalose cycle
            % sum Tre-cycle
            subplot(5,5,7)
            for j = 1:npSets
                if j == setup.refParams
                    sumTREcycle =  exp_G1P{j} + exp_T6P{j} + exp_TRE{j};
                    plot(dprofile{j}, sumTREcycle,'-o','color','black') 
                else
                    sumTREcycle =  exp_G1P{j} + exp_T6P{j} + exp_TRE{j};
                    plot(dprofile{j}, sumTREcycle, 'color', colorSet(j,:))
                end
            hold on
            end
            title('sum tre-cycle')
            % G1P 
            subplot(5,5,8)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_G1P{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_G1P{j}, 'color', colorSet(j,:))
                end
            hold on
            end
            title('G1P')
            % T6P
            subplot(5,5,9)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_T6P{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_T6P{j}, 'color', colorSet(j,:))
                end
            hold on
            end            
            title('T6P')
            % TRE
            subplot(5,5,10)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_TRE{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_TRE{j}, 'color', colorSet(j,:))
                end
            hold on
            end
            title('TRE')
        % AXP + Pi
            % sumAXP
            subplot(5,5,11)
            for j = 1:npSets
                if j == setup.refParams
                    sumAXP = exp_ATP{j} + exp_ADP{j} + exp_AMP{j};
                    plot(dprofile{j}, sumAXP,'-o','color','black') 
                else
                    sumAXP = exp_ATP{j} + exp_ADP{j} + exp_AMP{j};
                    plot(dprofile{j}, sumAXP, 'color', colorSet(j,:))
                end
            hold on
            end
            title('sumAXP')
            % ATP
            subplot(5,5,12)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_ATP{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_ATP{j}, 'color', colorSet(j,:))
                end
            hold on
            end
            title('ATP')
            % ADP
            subplot(5,5,13)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_ADP{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_ADP{j}, 'color', colorSet(j,:))
                end
            hold on
            end
            title('ADP')
            % AMP
            subplot(5,5,14)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_AMP{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_AMP{j}, 'color', colorSet(j,:))
                end
            hold on
            end
            title('AMP')
            % Pi
            subplot(5,5,15)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_PI{j},'-o','color','black') 
                else
                    plot(dprofile{j}, exp_PI{j}, 'color', colorSet(j,:))
                end
            hold on
            end
            title('Pi')
        % IXP
            % sumIXP
            subplot(5,5,17)
            for j = 1:npSets
                if j == setup.refParams
                    sumIXP = exp_IMP{j} + exp_INO{j} + exp_HYP{j};
                    plot(dprofile{j}, sumIXP,'-o','color','black') 
                else
                    sumIXP = exp_IMP{j} + exp_INO{j} + exp_HYP{j};
                    plot(dprofile{j}, sumIXP, 'color', colorSet(j,:))
                end
            hold on
            end
            title('umIXP')
            % IMP
            subplot(5,5,18)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_IMP{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_IMP{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('IMP')
            % INO
            subplot(5,5,19)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_INO{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_INO{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('INO')
            % HYP
            subplot(5,5,20)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_HYP{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_HYP{j},'color', colorSet(j,:))
                end
            hold on
            end            
            title('HYP')
        % NADX
            % sumNAD
            subplot(5,5,22)
            for j = 1:npSets
                if j == setup.refParams
                    sumNADX = exp_NAD{j} + exp_NADH{j};
                    plot(dprofile{j}, sumNADX,'-o','color','black') 
                else
                    sumNADX = exp_NAD{j} + exp_NADH{j};
                    plot(dprofile{j}, sumNADX, 'color', colorSet(j,:))
                end
            hold on
            end
            title('sumNADX')
            % NAD
            subplot(5,5,23)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_NAD{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_NAD{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('NAD')
            % NADH
            subplot(5,5,24)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_NADH{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_NADH{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('NADH')
        % set(gcf,'color','w'); suptitle('fig41. Metabolite pools. Steady state profile')
    end
    end
    
    % fig42. Metabolite pools. Fluxes related. Steady state profile
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Fluxes pools. Steady state profile','units','normalized','outerposition',[0 0 1 1])
        % C-structures
            % vGLT
            subplot(5,10,3)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_GLT{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_GLT{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vGLT')
            % vETtr
            subplot(5,10,4)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_ET_tr{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_ET_tr{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vEthanolTransport')
            % vGHtr
            subplot(5,10,5)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_GH_tr{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_GH_tr{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vGlycerolTransport')
            % vsinksTotal
            subplot(5,10,6)
            for j = 1:npSets
                if j == setup.refParams
                    exp_sinks_total = - exp_v_sinkPYR{j} - exp_v_sinkPEP{j} - exp_v_sinkP3G{j} + exp_v_sinkGAP{j} - exp_v_sinkG6P{j} + exp_v_sinkF6P{j} - exp_v_sinkACE{j};
                    plot(dprofile{j},exp_sinks_total','-o','color','black') 
                else
                    exp_sinks_total = - exp_v_sinkPYR{j} - exp_v_sinkPEP{j} - exp_v_sinkP3G{j} + exp_v_sinkGAP{j} - exp_v_sinkG6P{j} + exp_v_sinkF6P{j} - exp_v_sinkACE{j};
                    plot(dprofile{j},exp_sinks_total','color', colorSet(j,:))
                end
            hold on
            end
            title('vSinksTotal')
            % vtrecycle (pg1+tps1-nth1)
            subplot(5,10,2)
            for j = 1:npSets
                if j == setup.refParams
                    exp_tre_cycle = exp_v_PGM1{j} + exp_v_TPS1{j} - exp_v_NTH1{j};
                    plot(dprofile{j},exp_tre_cycle','-o','color','black') 
                else
                    exp_tre_cycle = exp_v_PGM1{j} + exp_v_TPS1{j} - exp_v_NTH1{j};
                    plot(dprofile{j},exp_tre_cycle','color', colorSet(j,:))
                end
            hold on
            end
            title('vTreCycle')
            % Balance C-fluxes
            subplot(5,10,1)
            for j = 1:npSets
                if j == setup.refParams
%                     vsumCstruct = exp_v_GLT{j} * 6 - exp_v_ET_tr{j} * 2 - exp_v_GH_tr{j} * 3 - (- exp_v_sinkPYR{j} * 3 - exp_v_sinkPEP{j} * 3 - exp_v_sinkP3G{j} * 3 + exp_v_sinkGAP{j} * 3 - exp_v_sinkG6P{j} * 6 + exp_v_sinkF6P{j} * 6 - exp_v_sinkACE{j} * 2) - exp_v_PDC{j};
                    mainflux = exp_v_GLT{j} * 6 - exp_v_G3PDH{j} * 3 - exp_v_ADH{j} * 2 - exp_v_PDC{j};
%                     sinksflux = -(exp_v_sinkG6P{j} * 6 + exp_v_sinkF6P{j} * 6 + exp_v_sinkGAP{j} * 3 + exp_v_sinkP3G{j} * 3 + exp_v_sinkPEP{j} * 3 + exp_v_sinkPYR{j} * 3 + exp_v_sinkACE{j} * 2);
                    sinksflux = -exp_v_sinkG6P{j} * 6 + exp_v_sinkF6P{j} * 6 + exp_v_sinkGAP{j} * 3 - exp_v_sinkP3G{j} * 3 - exp_v_sinkPEP{j} * 3 - exp_v_sinkPYR{j} * 3 - exp_v_sinkACE{j} * 2;
%                     vsumCstruct = mainflux - sinksflux;
                    vsumCstruct = mainflux + sinksflux;
                    plot(dprofile{j},vsumCstruct,'-o','color','black') 
                else
%                     vsumCstruct = exp_v_GLT{j} * 6 - exp_v_ET_tr{j} * 2 - exp_v_GH_tr{j} * 3 - (- exp_v_sinkPYR{j} * 3 - exp_v_sinkPEP{j} * 3 - exp_v_sinkP3G{j} * 3 + exp_v_sinkGAP{j} * 3 - exp_v_sinkG6P{j} * 6 + exp_v_sinkF6P{j} * 6 - exp_v_sinkACE{j} * 2) - exp_v_PDC{j};
                    mainflux = exp_v_GLT{j} * 6 - exp_v_G3PDH{j} * 3 - exp_v_ADH{j} * 2 - exp_v_PDC{j};
                    sinksflux = -(exp_v_sinkG6P{j} * 6 + exp_v_sinkF6P{j} * 6 + exp_v_sinkGAP{j} * 3 + exp_v_sinkP3G{j} * 3 + exp_v_sinkPEP{j} * 3 + exp_v_sinkPYR{j} * 3 + exp_v_sinkACE{j} * 2);
                    vsumCstruct = mainflux - sinksflux;
                    plot(dprofile{j},vsumCstruct,'color', colorSet(j,:))
                end
            hold on
            mainflux = canelas_SS.mtlD.v_GLT * 6 - canelas_SS.mtlD.v_G3PDH * 3 - canelas_SS.mtlD.v_ADH * 2 - canelas_SS.mtlD.v_PDC;
            sinksflux = canelas_SS.mtlD.sinkG6P * 6 + canelas_SS.mtlD.sinkF6P * 6 + canelas_SS.mtlD.sinkGAP * 3 + canelas_SS.mtlD.sinkP3G * 3 + canelas_SS.mtlD.sinkPEP * 3 + canelas_SS.mtlD.sinkPYR * 3 + canelas_SS.mtlD.sinkACE * 2;
            vsumCstruct_exp = mainflux - sinksflux;
            plot(dprofile{j},vsumCstruct_exp,'r+')
            end
            title('Balance C-fluxes')    % sum of c-struct fluxes
            % pgm1
            subplot(5,10,12)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_PGM1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_PGM1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vPGM1')
            % tps1
            subplot(5,10,13)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_TPS1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_TPS1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vTPS1')
            % tps2
            subplot(5,10,14)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_TPS2{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_TPS2{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vTPS2')
            % nth1
            subplot(5,10,15)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_NTH1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_NTH1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vNTH1')
        % AXP + Pi
            % ADK1
            subplot(5,10,21)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_ADK1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_ADK1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vADK1')
            % mito
            subplot(5,10,22)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_mito{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_mito{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vmito')
            % ATPase
            subplot(5,10,23)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_ATPase{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_ATPase{j},'color', colorSet(j,:))
                end
            hold on
            end            
            title('vATPase')
            % v_GLK
            subplot(5,10,24)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_HK{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_HK{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vGLK')
            % v_PFK
            subplot(5,10,25)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_PFK{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_PFK{j},'color', colorSet(j,:))
                end
            hold on
            end 
            title('vPFK')
            % v_PGK
            subplot(5,10,26)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_PGK{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_PGK{j},'color', colorSet(j,:))
                end
            hold on
            end 
            title('vPGK')
            % v_PYK
            subplot(5,10,27)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_PYK{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_PYK{j},'color', colorSet(j,:))
                end
            hold on
            end             
            title('vPYK')
            % v_TPS1
            subplot(5,10,28)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_TPS1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_TPS1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vTPS1')
            % v_Amd1
            subplot(5,10,29)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Amd1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Amd1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vAmd1')
            % v_Ade1312
            subplot(5,10,30)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Ade1213{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Ade1213{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vAde1312')
        % IXP
            % v_Hpt1
            subplot(5,10,32)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Hpt1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Hpt1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vHtp1')
            % v_Isn1
            subplot(5,10,33)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Isn1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Isn1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vIsn1')
            % v_Pnp1
            subplot(5,10,34)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Pnp1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Pnp1{j},'color', colorSet(j,:))
                end
            hold on
            end            
            title('vPnp1')
            % v_Amd1
            subplot(5,10,35)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Amd1{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Amd1{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vAmd1')
            % v_Ade1312
            subplot(5,10,36)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_Ade1213{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_Ade1213{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vAde1312')
        % NADX
            % v_G3PDH
            subplot(5,10,42)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_G3PDH{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_G3PDH{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vG3PDH')
            % v_GAPDH
            subplot(5,10,43)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_GAPDH{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_GAPDH{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vGAPDH')
            % v_ADH
            subplot(5,10,44)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_ADH{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_ADH{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vADH')
            % v_mitoNADH
            subplot(5,10,45)
            for j = 1:npSets
                if j == setup.refParams
                    plot(dprofile{j},exp_v_mitoNADH{j},'-o','color','black') 
                else
                    plot(dprofile{j},exp_v_mitoNADH{j},'color', colorSet(j,:))
                end
            hold on
            end
            title('vmitoNADH')
        % set(gcf,'color','w'); suptitle('fig42. Fluxes pools. Steady state profile')
    end
    end
    
    % fig43. Metabolite pools. Dynamic glucose pulse profile
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Metabolite pools. Dynamic glucose pulse profile','units','normalized','outerposition',[0 0 1 1])
        % C-structures
            % sum C-structures
            % glc,g6p,f6p,fbp,gap,dhap,g3p,glyc,g1p,t6p,tre,bpg,3pg,2pg,pep,pyr,ace,etoh
            subplot(5,5,2)
%             sumCstruct = Ygs(:,6) + Ygs(:,5) + Ygs(:,4) + Ygs(:,3) + Ygs(:,14) + Ygs(:,17) + Ygs(:,18) + Ygs(:,21) + Ygs(:,26) + Ygs(:,2) + Ygs(:,11) + Ygs(:,10) + Ygs(:,12) + Ygs(:,13) + Ygs(:,20) + Ygs(:,1) + Ygs(:,19) + Ygs(:,25);
%             plot(Tgs,sumCstruct,'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
%                     plot(Tgs{j}, Vgs{j}(:,2), 'k','LineWidth',1.2), 
                    sumCstruct = Ygs{j}(:,6) + Ygs{j}(:,5) + Ygs{j}(:,4) + Ygs{j}(:,3) + Ygs{j}(:,14) + Ygs{j}(:,17) + Ygs{j}(:,18) + Ygs{j}(:,21) + Ygs{j}(:,26) + Ygs{j}(:,2) + Ygs{j}(:,11) + Ygs{j}(:,10) + Ygs{j}(:,12) + Ygs{j}(:,13) + Ygs{j}(:,20) + Ygs{j}(:,1) + Ygs{j}(:,19) + Ygs{j}(:,25);
                    plot(Tgs{j},sumCstruct,'-','color','black') 
                    
                else
                    sumCstruct = Ygs{j}(:,6) + Ygs{j}(:,5) + Ygs{j}(:,4) + Ygs{j}(:,3) + Ygs{j}(:,14) + Ygs{j}(:,17) + Ygs{j}(:,18) + Ygs{j}(:,21) + Ygs{j}(:,26) + Ygs{j}(:,2) + Ygs{j}(:,11) + Ygs{j}(:,10) + Ygs{j}(:,12) + Ygs{j}(:,13) + Ygs{j}(:,20) + Ygs{j}(:,1) + Ygs{j}(:,19) + Ygs{j}(:,25);
                    plot(Tgs{j}, sumCstruct, 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('sum C-structures')
        % Trehalose cycle
            % sum Tre-cycle
            subplot(5,5,7)
%             sumTREcycle =  Ygs{j}(:,21) + Ygs{j}(:,26) + Ygs{j}(:,25);
%             plot(Tgs{j},sumTREcycle,'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
%                     plot(Tgs{j}, Vgs{j}(:,2), 'k','LineWidth',1.2), 
                    sumTREcycle =  Ygs{j}(:,21) + Ygs{j}(:,26) + Ygs{j}(:,25);
                    plot(Tgs{j},sumTREcycle,'-','color','black') 
                    
                else
                    sumTREcycle =  Ygs{j}(:,21) + Ygs{j}(:,26) + Ygs{j}(:,25);
                    plot(Tgs{j}, sumTREcycle, 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('sum tre-cycle')
            % G1P 
            subplot(5,5,8)
%             plot(Tgs{j},Ygs{j}(:,21),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,21),'-','color','black') 
                    
                else
                    plot(Tgs{j}, Ygs{j}(:,21), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('G1P')
            % T6P
            subplot(5,5,9)
%             plot(Tgs{j},Ygs{j}(:,26),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,26),'-','color','black') 
                    
                else
                    plot(Tgs{j},Ygs{j}(:,26), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('T6P')
            % TRE
            subplot(5,5,10)
%             plot(Tgs{j},Ygs{j}(:,25),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,25),'-','color','black') 
                    
                else
                    plot(Tgs{j},Ygs{j}(:,25), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('T6P')
            xlim([-100 340]), title('TRE')
        % AXP + Pi
            % sumAXP
            subplot(5,5,11)
%             sumAXP = Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16);
%             plot(Tgs{j},sumAXP,'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    sumAXP = Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16);
                    plot(Tgs{j},sumAXP,'-','color','black') 
                else
                    sumAXP = Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16);
                    plot(Tgs{j}, sumAXP, 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('sumAXP')
            % ATP
            subplot(5,5,12)
%             plot(Tgs{j},Ygs{j}(:,9),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,9),'-','color','black') 
                    
                else
                    plot(Tgs{j},Ygs{j}(:,9), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('ATP')
            % ADP
            subplot(5,5,13)
%             plot(Tgs{j},Ygs{j}(:,15),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,15),'-','color','black') 
                    
                else
                    plot(Tgs{j},Ygs{j}(:,15), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('ADP')
            % AMP
            subplot(5,5,14)
%             plot(Tgs{j},Ygs{j}(:,16),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,16),'-','color','black') 
                    
                else
                    plot(Tgs{j},Ygs{j}(:,16), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('AMP')
            % Pi
            subplot(5,5,15)
%             plot(Tgs{j},Ygs{j}(:,27),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,27),'-','color','black') 
                    
                else
                    plot(Tgs{j},Ygs{j}(:,27),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('Pi')
        % IXP
            % sumIXP
            subplot(5,5,17)
%             sumIXP = Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30);
%             plot(Tgs{j},sumIXP,'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    sumIXP = Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30);
                    plot(Tgs{j},sumIXP,'-','color','black') 
                    
                else
                    sumIXP = Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30);
                    plot(Tgs{j},sumIXP,'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('umIXP')
            % IMP
            subplot(5,5,18)
%             plot(Tgs{j},Ygs{j}(:,28),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,28),'-','color','black') 
                else
                    plot(Tgs{j},Ygs{j}(:,28),'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('IMP')
            % INO
            subplot(5,5,19)
%             plot(Tgs{j},Ygs{j}(:,29),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,29),'-','color','black') 
                else
                    plot(Tgs{j},Ygs{j}(:,29),'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('INO')
            % HYP
            subplot(5,5,20)
%             plot(Tgs{j},Ygs{j}(:,30),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,30),'-','color','black') 
                else
                    plot(Tgs{j},Ygs{j}(:,30),'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('HYP')
        % NADX
            % sumNAD
            subplot(5,5,22)
            sumNADX = Ygs{j}(:,7) + Ygs{j}(:,8);
            plot(Tgs{j},sumNADX,'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    sumNADX = Ygs{j}(:,7) + Ygs{j}(:,8);
                    plot(Tgs{j},sumNADX,'-','color','black') 
                else
                    sumNADX = Ygs{j}(:,7) + Ygs{j}(:,8);
                    plot(Tgs{j},sumNADX,'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('sumNADX')
            % NAD
            subplot(5,5,23)
%             plot(Tgs{j},Ygs{j}(:,7),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,7),'-','color','black') 
                else
                    plot(Tgs{j},Ygs{j}(:,7),'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('NAD')
            % NADH
            subplot(5,5,24)
%             plot(Tgs{j},Ygs{j}(:,8),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Ygs{j}(:,8),'-','color','black') 
                else
                    plot(Tgs{j},Ygs{j}(:,8),'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340]), title('NADH')
        % set(gcf,'color','w'); suptitle('fig43. Metabolite pools. Dynamic glucose pulse profile')
    end
    end
    
    % fig44. Metabolite pools. Fluxes related. Dynamic glucose pulse profile
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Metabolite pools. Fluxes related. Dynamic glucose pulse profile','units','normalized','outerposition',[0 0 1 1])
        % C-structures
            % vGLT
            subplot(5,10,2)
%             plot(Tgs{j},Vgs{j}(:,1),'-','color','black') 
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,1),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,1), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vGLT')
            % vETtr
            subplot(5,10,3)
%             plot(Tgs{j},Vgs{j}(:,23),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,23),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,23), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vEthanolTransport')
            % vGHtr
            subplot(5,10,4)
%             plot(Tgs{j},Vgs{j}(:,24),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,24),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,24),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vGlycerolTransport')
            % vsinksTotal
            subplot(5,10,5)
%             V_sinks_total = - Vgs{j}(:,39) - Vgs{j}(:,38) - Vgs{j}(:,37) + Vgs{j}(:,36) - Vgs{j}(:,34) + Vgs{j}(:,35) - Vgs{j}(:,40);
%             plot(Tgs{j},V_sinks_total,'-','color','black')
            for j = 1:npSets
                V_sinks_total = - Vgs{j}(:,39) - Vgs{j}(:,38) - Vgs{j}(:,37) + Vgs{j}(:,36) - Vgs{j}(:,34) + Vgs{j}(:,35) - Vgs{j}(:,40);
                if j == setup.refParams
                    plot(Tgs{j},V_sinks_total,'-','color','black') 
                    
                else
                    plot(Tgs{j},V_sinks_total,'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vSinksTotal')
            % vtrecycle (pg1+tps1-nth1)
            subplot(5,10,6)
%             V_tre_cycle = Vgs{j}(:,17) + Vgs{j}(:,21) - Vgs{j}(:,20);
%             plot(Tgs{j},V_tre_cycle,'-','color','black')
            for j = 1:npSets
                V_tre_cycle = Vgs{j}(:,17) + Vgs{j}(:,21) - Vgs{j}(:,20);
                if j == setup.refParams
                    plot(Tgs{j},V_tre_cycle,'-','color','black') 
                    
                else
                    plot(Tgs{j},V_tre_cycle,'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vTreCycle')
            % Balance C-fluxes
            subplot(5,10,7)
%             vsumCstruct = Vgs{j}(:,1) * 6 - Vgs{j}(:,23) * 2 - Vgs{j}(:,24) * 3 - (- 3 * Vgs{j}(:,39) - 3 * Vgs{j}(:,38) - 3 * Vgs{j}(:,37) + 3 * Vgs{j}(:,36) - 6 * Vgs{j}(:,34) + 6 * Vgs{j}(:,35) - 2 * Vgs{j}(:,40)) - Vgs{j}(:,13);
%             plot(Tgs{j},vsumCstruct,'-','color','black')
            for j = 1:npSets
                vsumCstruct = Vgs{j}(:,1) * 6 - Vgs{j}(:,23) * 2 - Vgs{j}(:,24) * 3 - (- 3 * Vgs{j}(:,39) - 3 * Vgs{j}(:,38) - 3 * Vgs{j}(:,37) + 3 * Vgs{j}(:,36) - 6 * Vgs{j}(:,34) + 6 * Vgs{j}(:,35) - 2 * Vgs{j}(:,40)) - Vgs{j}(:,13);
                if j == setup.refParams
                    plot(Tgs{j},vsumCstruct,'-','color','black') 
                    
                else
                    plot(Tgs{j},vsumCstruct,'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('Balance C-fluxes \color{red}(multiply each flux by the carbon structure)')    % sum of c-struct fluxes
            % pgm1
            subplot(5,10,12)
%             plot(Tgs{j},Vgs{j}(:,17),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,17),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,17),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vPGM1')
            % tps1
            subplot(5,10,13)
%             plot(Tgs{j},Vgs{j}(:,21),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,21),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,21),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vTPS1')
            % tps2
            subplot(5,10,14)
%             plot(Tgs{j},Vgs{j}(:,19),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,19),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,19),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vTPS2')
            % nth1
            subplot(5,10,15)
%             plot(Tgs{j},Vgs{j}(:,20),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,20),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,20),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vNTH1')
        % AXP + Pi
            % ADK1
            subplot(5,10,21)
%             plot(Tgs{j},Vgs{j}(:,14),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,14),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,14),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vADK1')
            % mito
            subplot(5,10,22)
%             plot(Tgs{j},Vgs{j}(:,28),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,28),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,28),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vmito')
            % ATPase
            subplot(5,10,23)
%             plot(Tgs{j},Vgs{j}(:,27),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,27),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,27),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vATPase')
            % v_GLK
            subplot(5,10,24)
%             plot(Tgs{j},Vgs{j}(:,2),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,2),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,2),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vGLK')
            % v_PFK
            subplot(5,10,25)
%             plot(Tgs{j},Vgs{j}(:,4),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,4),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,4),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vPFK')
            % v_PGK
            subplot(5,10,26)
%             plot(Tgs{j},Vgs{j}(:,9),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,9),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,9),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vPGK')
            % v_PYK
            subplot(5,10,27)
%             plot(Tgs{j},Vgs{j}(:,12),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,12),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,12),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vPYK')
            % v_TPS1
            subplot(5,10,28)
%             plot(Tgs{j},Vgs{j}(:,21),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,21),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,21),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vTPS1')
            % v_Amd1
            subplot(5,10,29)
%             plot(Tgs{j},Vgs{j}(:,29),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,29),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,29),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vAmd1')
            % v_Ade1312
            subplot(5,10,30)
%             plot(Tgs{j},Vgs{j}(:,30),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,30),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,30),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vAde1312')
        % IXP
            % v_Hpt1
            subplot(5,10,32)
%             plot(Tgs{j},Vgs{j}(:,33),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,33),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,33),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vHtp1')
            % v_Isn1
            subplot(5,10,33)
%             plot(Tgs{j},Vgs{j}(:,31),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,31),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,31),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vIsn1')
            % v_Pnp1
            subplot(5,10,34)
%             plot(Tgs{j},Vgs{j}(:,32),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,32),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,32),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vPnp1')
            % v_Amd1
            subplot(5,10,35)
%             plot(Tgs{j},Vgs{j}(:,29),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,29),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,29),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vAmd1')
            % v_Ade1312
            subplot(5,10,36)
%             plot(Tgs{j},Vgs{j}(:,30),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,30),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,30),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vAde1312')
        % NADX
            % v_G3PDH
            subplot(5,10,42)
%             plot(Tgs{j},Vgs{j}(:,6),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,6),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,6),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vG3PDH')
            % v_GAPDH
            subplot(5,10,43)
%             plot(Tgs{j},Vgs{j}(:,8),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,8),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,8),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vGAPDH')
            % v_ADH
            subplot(5,10,44)
%             plot(Tgs{j},Vgs{j}(:,41),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,41),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,41),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vADH')
            % v_mitoNADH
            subplot(5,10,45)
%             plot(Tgs{j},Vgs{j}(:,26),'-','color','black')
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j},Vgs{j}(:,26),'-','color','black') 
                    
                else
                    plot(Tgs{j},Vgs{j}(:,26),'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340]), title('vmitoNADH')
        % set(gcf,'color','w'); suptitle('fig44. Metabolite pools. Fluxes related. Dynamic glucose pulse profile')
    end
    end

end

% old style plots profiles
if(setup.plotResultsMode == 40)

    % fig1. Steady state metabolite profile. All
    for i = 1
    if setup.runSScanelas == 1

    figure('name','Steady state metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
    % figure(101)
    % G6P
    subplot(4,7,1)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.G6P, 'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_G6P{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_G6P{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('G6P')
    xlabel('growth rate [h-1]')
    ylabel('concentration [mM]')
    % xlim([0 0.2])

    % F6P
    subplot(4,7,2)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.F6P,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_F6P{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_F6P{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('F6P')
    % xlim([0 0.2])

    % FBP
    subplot(4,7,3)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.FBP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_FBP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_FBP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('FBP')
    % xlim([0 0.2])

    % GAP
    subplot(4,7,4)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.GAP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_GAP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_GAP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('GAP')
    % xlim([0 0.2])

    % DHAP
    subplot(4,7,5)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.DHAP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_DHAP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_DHAP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('DHAP')
    % xlim([0 0.2])

    % G1P
    subplot(4,7,6)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.G1P,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_G1P{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_G1P{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('G1P')
    % xlim([0 0.2])

    % T6P
    subplot(4,7,7)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.T6P,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_T6P{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_T6P{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('T6P')
    % xlim([0 0.2])

    % G3P
    subplot(4,7,8)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.G3P,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_G3P{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_G3P{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('G3P')
    % xlim([0 0.2])

    % P3G
    subplot(4,7,9)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.P3G,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_P3G{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_P3G{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('P3G')
    % xlim([0 0.2])

    % P2G
    subplot(4,7,10)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.P2G,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_P2G{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_P2G{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('P2G')
    % xlim([0 0.2])

    % PEP
    subplot(4,7,11)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.PEP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_PEP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_PEP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('PEP')
    % xlim([0 0.2])

    % PYR
    subplot(4,7,12)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.PYR,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_PYR{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_PYR{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('PYR')
    % xlim([0 0.2])

    % ETOH
    subplot(4,7,13)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.ETOH,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_ETOH{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_ETOH{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('ETOH')
    % xlim([0 0.2])

    % AMP
    subplot(4,7,14)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.AMP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_AMP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_AMP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('AMP')
    % xlim([0 0.2])

    % ADP
    subplot(4,7,15)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.ADP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_ADP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_ADP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('ADP')
    % xlim([0 0.2])

    % ATP
    subplot(4,7,16)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.ATP,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_ATP{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_ATP{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('ATP')
    % xlim([0 0.2])

    % NADH
    subplot(4,7,17)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.NADH,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_NADH{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_NADH{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('NADH')
    % xlim([0 0.2])

    % NAD
    subplot(4,7,18)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.NAD,'r+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_NAD{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_NAD{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('NAD')
    % xlim([0 0.2])
    Qo = legend('Simulation','Experimental data');
    set(Qo, 'Position', [0.025,0.85,0.075,0.10]);

        % Latest additions
    % Pi
    subplot(4,7,19)
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_PI{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_PI{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('Pi')
    % xlim([0 0.2])

    % BPG
    subplot(4,7,20)
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_BPG{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_BPG{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('BPG')
    % xlim([0 0.2])

    % GLCi
    subplot(4,7,21)
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_GLCi{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_GLCi{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('GLCi')
    % xlim([0 0.2])

    % ACE
    subplot(4,7,22)
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_ACE{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_ACE{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('ACE')
    % xlim([0 0.2])
    
    subplot(4,7,23)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.TRE, 'b+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_TRE{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_TRE{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('TRE')
    xlabel('growth rate [h-1]')
    ylabel('concentration [mM]')
    % xlim([0 0.2])
    
    subplot(4,7,24)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.ECetoh, 'b+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_ECetoh{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_ECetoh{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('Etoh_{EC}')
    xlabel('growth rate [h-1]')
    ylabel('concentration [mM]')
    % xlim([0 0.2])
    
    subplot(4,7,25)
    plot(canelas_SS.mtlD.D, canelas_SS.mtlD.ECglycerol, 'b+')
    hold on
    for j = 1:npSets
        if j == setup.refParams
            plot(dprofile{j}(1:q(1)), exp_ECglycerol{j}, 'k','LineWidth',1.2)
        else
            plot(dprofile{j}(1:q(1)), exp_ECglycerol{j}, 'color', colorSet(j,:))
        end
        hold on
    end
    title('Glyc_{EC}')
    xlabel('growth rate [h-1]')
    ylabel('concentration [mM]')
    % xlim([0 0.2])

    % set(gcf,'color','w'); suptitle('fig1. Steady state metabolite profile. All')
    end
    end

    % fig2. Steady state flux profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state flux profile. All','units','normalized','outerposition',[0 0 1 1])

        subplot(4,8,1)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GLT,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('GLT')
        xlabel('growth rate [h-1]')
        ylabel('reaction rate [mM s-1]')
        % xlim([0 0.2])

        % v_HK
        subplot(4,8,2)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_HK,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_HK{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_HK{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('HK')
        % xlim([0 0.2])

        % v_PGI
        subplot(4,8,3)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PGI,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PGI{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PGI{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        % ylim([0 2]), 
        title('PGI')
        % xlim([0 0.2])

        % v_PFK
        subplot(4,8,4)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PFK, 'r+', 'MarkerSize', 5)
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PFK{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PFK{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('PFK')
        % xlim([0 0.2])

        % v_FBA
        subplot(4,8,5)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_FBA,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_FBA{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_FBA{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('FBA')
        % xlim([0 0.2])

        % v_TPI
        subplot(4,8,6)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_TPI,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_TPI{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_TPI{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('TPI')
        % xlim([0 0.2])

        % v_G3PDH
        subplot(4,8,7)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_G3PDH,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_G3PDH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_G3PDH{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('G3PDH')
        % xlim([0 0.2])

        % v_GAPDH
        subplot(4,8,8)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GAPDH,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GAPDH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GAPDH{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('GAPDH')
        % xlim([0 0.2])

        % v_PGK
        subplot(4,8,9)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PGK,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PGK{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PGK{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGK')
        % xlim([0 0.2])

        % v_PGM
        subplot(4,8,10)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PGM,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PGM{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PGM{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGM')
        % xlim([0 0.2])

        % v_ENO
        subplot(4,8,11)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ENO,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ENO{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ENO{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('ENO')
        % xlim([0 0.2])

        % v_PYK
        subplot(4,8,12)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PYK,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PYK{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PYK{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('PYK')
        % xlim([0 0.2])

        % v_PDC
        subplot(4,8,13)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_PDC,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PDC{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PDC{j}, 'color', colorSet(j,:))
            end
            hold on
        end
%         hold on
        title('PDC')
        % xlim([0 0.2])

        % v_ADH
        subplot(4,8,14)
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH,'r+')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ADH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ADH{j}, 'color', colorSet(j,:))
            end
            hold on
        end
%         hold on
        title('ADH')
        % xlim([0 0.2])

        subplot(4,8,15)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PDH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PDH{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('PDH')

        subplot(4,8,16)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ACE{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ACE{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('ACE')

        subplot(4,8,17)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.qEtoh,'b+')
        title('ET_tr')

        subplot(4,8,18)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_HOR2{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_HOR2{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('HOR2')

        subplot(4,8,19)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'color', colorSet(j,:))
            end
            hold on
            plot(canelas_SS.mtlD.D, canelas_SS.mtlD.qGlyc,'b+')
        end
        title('GH_tr')

        subplot(4,8,20)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_PGM1{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_PGM1{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGM1')

        subplot(4,8,21)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_TPS1{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_TPS1{j}, 'color', colorSet(j,:))
            end
            hold on
        end        
        title('TPS1')

        subplot(4,8,22)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_TPS2{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_TPS2{j}, 'color', colorSet(j,:))
            end
            hold on
        end           
        title('TPS2')

        subplot(4,8,23)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_NTH1{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_NTH1{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('NTH1')

        d = dprofile{j}(1:q(1));
        poly_sinkG6P    = 4.2647 * d.^3 -   1.7552 * d.^2 -  0.5809 * d    - 0.0047;
        poly_sinkF6P    = 230.0029 * d.^6 - 252.2129 * d.^5 + 92.6247 * d.^4 - 14.1235 * d.^3 + 1.0306 * d.^2 + 0.2117 * d - 0.0005;
        poly_sinkGAP    = 184.4055 * d.^6 - 199.5226 * d.^5 + 73.2004 * d.^4 - 11.1086 * d.^3 + 0.7643 * d.^2 + 0.1234 * d - 0.00001;
        poly_sinkP3G    = -   0.0877 * d.^2 -   0.0892 * d   +  0.0004;
        poly_sinkPEP    = -   0.0627 * d.^2 -   0.0599 * d   +  0.0001;
        poly_sinkPYR    = - 7.4563e+03 * d.^6 + 8.2025e+03 * d.^5 - 3.2635e+03 * d.^4 + 581.4912 * d.^3 - 46.6358 * d.^2 - 0.0384 * d - 0.0226;
        poly_sinkACE    =     - 116.4277 * d.^6 -  20.8956 * d.^5 + 62.6872 * d.^4 - 24.9573 * d.^3 + 3.9700 * d.^2 - 0.5698 * d + 0.0039;

        subplot(4,8,24)
        plot(d, poly_sinkG6P, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkG6P{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkG6P{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('sinkG6P')

        subplot(4,8,25)
        plot(d, poly_sinkF6P, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkF6P{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkF6P{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('sinkF6P')

        subplot(4,8,26)
        plot(d, poly_sinkGAP, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkGAP{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkGAP{j}, 'color', colorSet(j,:))
            end
            hold on
        end
        title('sinkGAP')

        subplot(4,8,27)
        plot(d, poly_sinkP3G, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkP3G{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkP3G{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('sinkP3G')

        subplot(4,8,28)
        plot(d, poly_sinkPEP, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkPEP{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkPEP{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('sinkPEP')

        subplot(4,8,29)
        plot(d, poly_sinkPYR, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('sinkPYR')

        subplot(4,8,30)
        plot(d, poly_sinkACE, 'r-')
        hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_sinkACE{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_sinkACE{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('sinkACE')

        subplot(4,8,31)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_mitoNADH{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_mitoNADH{j}, 'color', colorSet(j,:))
            end
            hold on
        end    
        title('v_{mitoNADH}')
        
        subplot(4,8,32)
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_vacPi{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_vacPi{j}, 'color', colorSet(j,:))
            end
            hold on
        end           
        title('vac_{Pi}')      
                
        % set(gcf,'color','w'); suptitle('fig2. Steady state flux profile. All')
    end
    end

    %fig3. Dynamic glucose pulse metabolite profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
%         for k = 1:length(legenda.metabolites)
%             subplot(4,8,k)
%             % plot simulations
%             for j = 1:npSets
%                 if j == setup.refParams
%                     plot(Tgs{j}, Ygs{j}(:,k), 'k')
%                 else
%                     plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
%                 end
%                 hold on
%             end
%             % plot experimental data
%             if data.ordered.metabolites.onoff(k) == 1
%                 plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
%             end
%             title(legenda.metabolites(k))
%             xlim([-100 340])
%         end
        
        subplot(4,6,1)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,5), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,5), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,2),'r+')
        title('G6P')
        xlabel('time [s]')
        ylabel('concentration [mM]')
        xlim([-100 340])

        % F6P
        subplot(4,6,2)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,4), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,4), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,3),'r+')
        title('F6P')
        xlim([-100 340])

        % FBP
        subplot(4,6,3)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,3), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,3), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,4),'r+')
        title('FBP')
        xlim([-100 340])

        % GAP
        subplot(4,6,4)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,14), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,14), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,5),'r+')
        title('GAP')
        xlim([-100 340])

        % DHAP
        subplot(4,6,5)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,17), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,17), 'color', colorSet(j,:))
            end
            hold on
        end
        title('DHAP')
        xlim([-100 340])

        % G1P
        subplot(4,6,6)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,21), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,21), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,13),'r+')
        title('G1P')
        xlim([-100 340])

        % T6P
        subplot(4,6,7)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,26), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,26), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,15),'r+')
        title('T6P')
        xlim([-100 340])

        % G3P
        subplot(4,6,8)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,18), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,18), 'color', colorSet(j,:))
            end
            hold on
        end
        title('G3P')
        xlim([-100 340])

        % P3G
        subplot(4,6,9)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,11), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,11), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,6),'r+')
        title('P3G')
        xlim([-100 340])

        % P2G
        subplot(4,6,10)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,10), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,10), 'color', colorSet(j,:))
            end
            hold on
        end
        title('P2G')
        xlim([-100 340])

        % PEP
        subplot(4,6,11)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,12), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,12), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,7),'r+')
        title('PEP')
        xlim([-100 340])

        % PYR
        subplot(4,6,12)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,13), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,13), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PYR')
        xlim([-100 340])

        % ETOH
        subplot(4,6,13)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,20), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,20), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ETOH')
        xlim([-100 340])

        % AMP
        subplot(4,6,14)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,16), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,16), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_nucleotides,data.nucleotides(:,2),'r+')
        title('AMP')
        xlim([-100 340])

        % ADP
        subplot(4,6,15)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,15), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,15), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_nucleotides,data.nucleotides(:,3),'r+')
        title('ADP')
        xlim([-100 340])

        % ATP
        subplot(4,6,16)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,9), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,9), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_nucleotides,data.nucleotides(:,4),'r+')
        title('ATP')

        xlim([-100 340])
        Qo = legend('Simulation','Experimental data');
        set(Qo, 'Position', [0.025,0.85,0.075,0.10]);

        % NADH
        subplot(4,6,17)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,8), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,8), 'color', colorSet(j,:))
            end
            hold on
        end
        title('NADH')
        xlim([-100 340])

        % NAD
        subplot(4,6,18)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,7), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,7), 'color', colorSet(j,:))
            end
            hold on
        end
        title('NAD')
        xlim([-100 340])

        % Pi
        subplot(4,6,19)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,27), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,27), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Pi')
        xlim([-100 340])

        % BPG
        subplot(4,6,20)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,2), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,2), 'color', colorSet(j,:))
            end
            hold on
        end
        title('BPG')
        xlim([-100 340])

        % GLCi
        subplot(4,6,21)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,6), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,6), 'color', colorSet(j,:))
            end
            hold on
        end
        title('GLCi')
        xlim([-100 340])

        % ACE
        subplot(4,6,22)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,1), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,1), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ACE')
        xlim([-100 340])

        % TRE
        subplot(4,6,23)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,25), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,25), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites, data.metabolites(:,16),'r+')
        title('TRE')
        xlim([-100 340])

        % UDP_glc
        subplot(4,6,24)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,24), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,24), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites, data.metabolites(:,14),'r+')
        title('UDP_glc')
        xlim([-100 340])
        
        % set(gcf,'color','w'); suptitle('fig3. Dynamic glucose pulse metabolite profile. All')
    end
    end
    
    %fig4. Dynamic glucose pulse flux profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        subplot(4,8,1)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,1), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,1), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,1), 'r+')
        title('GLT')
        xlabel('time [s^{-1}]')
        ylabel('reaction rate [mM s^{-1}]')
        xlim([-100 340])

        % v_HK
        subplot(4,8,2)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,2), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,2), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,1), 'r+')
        title('HK')
        xlim([-100 340])

        % v_PGI
        subplot(4,8,3)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,3), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,3), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,2)-data.fluxes(:,3), 'r+')
        hold on
        title('PGI')
        xlim([-100 340])

        % v_PFK
        subplot(4,8,4)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,4), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,4), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,4), 'r+')
        hold on
        title('PFK')
        xlim([-100 340])

        % v_FBA
        subplot(4,8,5)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,5), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,5), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,5), 'r+')
        hold on
        title('FBA')
        xlim([-100 340])

        % v_TPI
        subplot(4,8,6)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,7), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,7), 'color', colorSet(j,:))
            end
            hold on
        end
        title('TPI')
        xlim([-100 340])

        % v_G3PDH
        subplot(4,8,7)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,6), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,6), 'color', colorSet(j,:))
            end
            hold on
        end
        title('G3PDH')
        xlim([-100 340])

        % v_GAPDH
        subplot(4,8,8)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,8), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,8), 'color', colorSet(j,:))
            end
            hold on
        end
        title('GAPDH')
        xlim([-100 340])

        % v_PGK
        subplot(4,8,9)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,9), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,9), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGK')
        xlim([-100 340])

        % v_PGM
        subplot(4,8,10)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,10), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,10), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGM')
        xlim([-100 340])

        % v_ENO
        subplot(4,8,11)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,11), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,11), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ENO')
        xlim([-100 340])

        % v_PYK
        subplot(4,8,12)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,12), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,12), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PYK')
        xlim([-100 340])

        % v_PDC
        subplot(4,8,13)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,13), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,13), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PDC')
        xlim([-100 340])

        % v_ADH
        subplot(4,8,14)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,41), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,41), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ADH')
        xlim([-100 340])

        % v_PDH
        subplot(4,8,15)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,25), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,25), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PDH')
        xlim([-100 340])

        % v_ACE
        subplot(4,8,16)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,22), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,22), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ACE')
        xlim([-100 340])

        % v_ET_tr
        subplot(4,8,17)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,23), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,23), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ET_tr')
        xlim([-100 340])

        % v_HOR2
        subplot(4,8,18)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,15), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,15), 'color', colorSet(j,:))
            end
            hold on
        end
        title('HOR2')
        xlim([-100 340])

        % v_GL_tr
        subplot(4,8,19)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,24), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,24), 'color', colorSet(j,:))
            end
            hold on
        end
        title('GL_tr')
        xlim([-100 340])

        % v_PGM1
        subplot(4,8,20)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, -Vgs{j}(:,17), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, -Vgs{j}(:,17), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,6)-data.fluxes(:,7), 'r+')
        title('PGM1')
        xlim([-100 340])

        % v_TPS1
        subplot(4,8,21)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,21), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,21), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,9), 'r+')
        title('TPS1')
        xlim([-100 340])

        % v_TPS2
        subplot(4,8,22)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,19), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,19), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,10), 'r+')
        title('TPS2')
        xlim([-100 340])

        % v_NTH1
        subplot(4,8,23)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,20), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,20), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,11), 'r+')
        title('NTH1')
        xlim([-100 340])

        % v_sinkG6P
        subplot(4,8,24)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,34), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,34), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{G6P}')
        xlim([-100 340])

        % v_sinkF6P
        subplot(4,8,25)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,35), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,35), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{F6P}')
        xlim([-100 340])

        % v_sinkGAP
        subplot(4,8,26)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,36), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,36), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{GAP}')
        xlim([-100 340])

        % v_sinkP3G
        subplot(4,8,27)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,37), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,37), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{P3G}')
        xlim([-100 340])

        % v_sinkPEP
        subplot(4,8,28)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,38), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,38), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{PEP}')
        xlim([-100 340])

        % v_sinkPYR
        subplot(4,8,29)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,39), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,39), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{PYR}')
        xlim([-100 340])

        % v_sinkACE
        subplot(4,8,30)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,40), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,40), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{ACE}')
        xlim([-100 340])

        %recently added (v_mitoNAD)
        subplot(4,8,31)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,26), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,26), 'color', colorSet(j,:))
            end
            hold on
        end
        title('v_{mitoNAD}')
        xlim([-100 340])
        
        subplot(4,8,32)
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,42), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,42), 'color', colorSet(j,:))
            end
            hold on
        end           
        title('vac_{Pi}')   
        xlim([-100 340])
        
        % set(gcf,'color','w'); suptitle('fig4. Dynamic glucose pulse flux profile. All')
    end
    end
    
    %fig5. Dynamic glucose pulse metabolite profile. All. Until 2000s
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse metabolite profile. All. Until 2000s','units','normalized','outerposition',[0 0 1 1])
        subplot(4,6,1)
%         plot(Tgs, Ygs(:,5), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,5), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,5), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,2),'r+')
        hold on
        title('G6P')
        xlabel('time [s]')
        ylabel('concentration [mM]')
        xlim([-100 2000])

        % F6P
        subplot(4,6,2)
%         plot(Tgs, Ygs(:,4), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,4), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,4), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,3),'r+')
        hold on
        title('F6P')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % FBP
        subplot(4,6,3)
%         plot(Tgs, Ygs(:,3), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,3), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,3), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,4),'r+')
        hold on
        title('FBP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % GAP
        subplot(4,6,4)
%         plot(Tgs, Ygs(:,14), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,14), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,14), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,5),'r+')
        hold on
        title('GAP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % DHAP
        subplot(4,6,5)
%         plot(Tgs, Ygs(:,17), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,17), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,17), 'color', colorSet(j,:))
            end
            hold on
        end
        title('DHAP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % G1P
        subplot(4,6,6)
%         plot(Tgs, Ygs(:,21), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,21), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,21), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,13),'r+')
        hold on
        title('G1P')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % T6P
        subplot(4,6,7)
%         plot(Tgs, Ygs(:,26), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,26), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,26), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,15),'r+')
        hold on
        title('T6P')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % G3P
        subplot(4,6,8)
%         plot(Tgs, Ygs(:,18), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,18), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,18), 'color', colorSet(j,:))
            end
            hold on
        end
        title('G3P')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % P3G
        subplot(4,6,9)
%         plot(Tgs, Ygs(:,11), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,11), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,11), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,6),'r+')
        hold on
        title('P3G')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % P2G
        subplot(4,6,10)
%         plot(Tgs, Ygs(:,10), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,10), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,10), 'color', colorSet(j,:))
            end
            hold on
        end
        title('P2G')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % PEP
        subplot(4,6,11)
%         plot(Tgs, Ygs(:,12), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,12), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,12), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites,data.metabolites(:,7),'r+')
        hold on
        title('PEP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % PYR
        subplot(4,6,12)
%         plot(Tgs, Ygs(:,13), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,13), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,13), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PYR')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % ETOH
        subplot(4,6,13)
%         plot(Tgs, Ygs(:,20), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,20), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,20), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ETOH')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % AMP
        subplot(4,6,14)
%         plot(Tgs, Ygs(:,16), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,16), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,16), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_nucleotides,data.nucleotides(:,2),'r+')
        hold on
        title('AMP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % ADP
        subplot(4,6,15)
%         plot(Tgs, Ygs(:,15), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,15), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,15), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_nucleotides,data.nucleotides(:,3),'r+')
        hold on
        title('ADP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % ATP
        subplot(4,6,16)
%         plot(Tgs, Ygs(:,9), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,9), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,9), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_nucleotides,data.nucleotides(:,4),'r+')
        hold on
        title('ATP')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])
        Qo = legend('Simulation','Experimental data');
        set(Qo, 'Position', [0.025,0.85,0.075,0.10]);

        % NADH
        subplot(4,6,17)
%         plot(Tgs, Ygs(:,8), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,8), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,8), 'color', colorSet(j,:))
            end
            hold on
        end
        title('NADH')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % NAD
        subplot(4,6,18)
%         plot(Tgs, Ygs(:,7), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,7), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,7), 'color', colorSet(j,:))
            end
            hold on
        end
        title('NAD')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

            % Latest additions

        % Pi
        subplot(4,6,19)
%         plot(Tgs, Ygs(:,27), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,27), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,27), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Pi')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % BPG
        subplot(4,6,20)
%         plot(Tgs, Ygs(:,2), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,2), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,2), 'color', colorSet(j,:))
            end
            hold on
        end
        title('BPG')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % GLCi
        subplot(4,6,21)
%         plot(Tgs, Ygs(:,6), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,6), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,6), 'color', colorSet(j,:))
            end
            hold on
        end
        title('GLCi')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % ACE
        subplot(4,6,22)
%         plot(Tgs, Ygs(:,1), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,1), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,1), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ACE')
        % xlabel('time [s]')
        % ylabel('concentration [mM]')
        xlim([-100 2000])

        % TRE
        subplot(4,6,23)
%         plot(Tgs, Ygs(:,25), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,25), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,25), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites, data.metabolites(:,16),'r+')
        hold on
        title('TRE')
        xlim([-100 2000])

        % UDP_glc
        subplot(4,6,24)
%         plot(Tgs, Ygs(:,24), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Ygs{j}(:,24), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Ygs{j}(:,24), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_metabolites, data.metabolites(:,14),'r+')
        hold on
        title('UDP_glc')
        xlim([-100 2000])
        % set(gcf,'color','w'); suptitle('fig5. Dynamic glucose pulse metabolite profile. All. Until 2000s')
    end
    end
    
    %fig6. Dynamic glucose pulse flux profile. All. Until 2000s
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        subplot(4,8,1)
%         plot(Tgs, Vgs(:,1), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,1), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,1), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,1), 'r+')
        hold on
        title('GLT')
        xlabel('time [s^{-1}]')
        ylabel('reaction rate [mM s^{-1}]')
        xlim([-100 2000])

        % v_HK
        subplot(4,8,2)
%         plot(Tgs, Vgs(:,2), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,2), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,2), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,1), 'r+')
        hold on
        title('HK')
        xlim([-100 2000])

        % v_PGI
        subplot(4,8,3)
%         plot(Tgs, Vgs(:,3), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,3), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,3), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,2)-data.fluxes(:,3), 'r+')
        hold on
        title('PGI')
        xlim([-100 2000])

        % v_PFK
        subplot(4,8,4)
%         plot(Tgs, Vgs(:,4), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,4), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,4), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,4), 'r+')
        hold on
        title('PFK')
        xlim([-100 2000])

        % v_FBA
        subplot(4,8,5)
%         plot(Tgs, Vgs(:,5), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,5), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,5), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(1:end),data.fluxes(:,5), 'r+')
        hold on
        title('FBA')
        xlim([-100 2000])

        % v_TPI
        subplot(4,8,6)
%         plot(Tgs, Vgs(:,7), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,7), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,7), 'color', colorSet(j,:))
            end
            hold on
        end
        title('TPI')
        xlim([-100 2000])

        % v_G3PDH
        subplot(4,8,7)
%         plot(Tgs, Vgs(:,6), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,6), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,6), 'color', colorSet(j,:))
            end
            hold on
        end
        title('G3PDH')
        xlim([-100 2000])

        % v_GAPDH
        subplot(4,8,8)
%         plot(Tgs, Vgs(:,8), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,8), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,8), 'color', colorSet(j,:))
            end
            hold on
        end
        title('GAPDH')
        xlim([-100 2000])

        % v_PGK
        subplot(4,8,9)
%         plot(Tgs, Vgs(:,9), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,9), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,9), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGK')
        xlim([-100 2000])

        % v_PGM
        subplot(4,8,10)
%         plot(Tgs, Vgs(:,10), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,10), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,10), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PGM')
        xlim([-100 2000])

        % v_ENO
        subplot(4,8,11)
%         plot(Tgs, Vgs(:,11), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,11), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,11), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ENO')
        xlim([-100 2000])

        % v_PYK
        subplot(4,8,12)
%         plot(Tgs, Vgs(:,12), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,12), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,12), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PYK')
        xlim([-100 2000])

        % v_PDC
        subplot(4,8,13)
%         plot(Tgs, Vgs(:,13), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,13), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,13), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PDC')
        xlim([-100 2000])

        % v_ADH
        subplot(4,8,14)
%         plot(Tgs, Vgs(:,41), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,41), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,41), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ADH')
        xlim([-100 2000])

        % v_PDH
        subplot(4,8,15)
%         plot(Tgs, Vgs(:,25), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,25), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,25), 'color', colorSet(j,:))
            end
            hold on
        end
        title('PDH')
        xlim([-100 2000])

        % v_ACE
        subplot(4,8,16)
%         plot(Tgs, Vgs(:,22), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,22), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,22), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ACE')
        xlim([-100 2000])

        % v_ET_tr
        subplot(4,8,17)
%         plot(Tgs, Vgs(:,23), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,23), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,23), 'color', colorSet(j,:))
            end
            hold on
        end
        title('ET_tr')
        xlim([-100 2000])

        % v_HOR2
        subplot(4,8,18)
%         plot(Tgs, Vgs(:,15), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,15), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,15), 'color', colorSet(j,:))
            end
            hold on
        end
        title('HOR2')
        xlim([-100 2000])

        % v_GL_tr
        subplot(4,8,19)
%         plot(Tgs, Vgs(:,24), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,24), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,24), 'color', colorSet(j,:))
            end
            hold on
        end
        title('GL_tr')
        xlim([-100 2000])

        % v_PGM1
        subplot(4,8,20)
%         if setup.clamp.TRE == 1
%             plot(Tgs, Vgs(:,17), 'k')
%         elseif setup.clamp.TRE == 0
%             plot(Tgs, -Vgs(:,17), 'k')
%         end
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, -Vgs{j}(:,17), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, -Vgs{j}(:,17), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,6)-data.fluxes(:,7), 'r+')
        hold on
        title('PGM1')
        xlim([-100 2000])

        % v_TPS1
        subplot(4,8,21)
%         plot(Tgs, Vgs(:,21), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,21), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,21), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,9), 'r+')
        hold on
        title('TPS1')
        xlim([-100 2000])

        % v_TPS2
        subplot(4,8,22)
%         plot(Tgs, Vgs(:,19), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,19), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,19), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,10), 'r+')
        hold on
        title('TPS2')
        xlim([-100 2000])

        % v_NTH1
        subplot(4,8,23)
%         plot(Tgs, Vgs(:,20), 'k')
%         hold on
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,20), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,20), 'color', colorSet(j,:))
            end
            hold on
        end
        plot(data.time_fluxes(:), data.fluxes(:,11), 'r+')
        hold on
        title('NTH1')
        xlim([-100 2000])

       % v_sinkG6P
        subplot(4,8,24)
%         plot(Tgs, Vgs(:,34), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,34), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,34), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{G6P}')
        xlim([-100 2000])

        % v_sinkF6P
        subplot(4,8,25)
%         plot(Tgs, Vgs(:,35), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,35), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,35), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{F6P}')
        xlim([-100 2000])

        % v_sinkGAP
        subplot(4,8,26)
%         plot(Tgs, Vgs(:,36), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,36), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,36), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{GAP}')
        xlim([-100 2000])

        % v_sinkP3G
        subplot(4,8,27)
%         plot(Tgs, Vgs(:,37), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,37), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,37), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{P3G}')
        xlim([-100 2000])

        % v_sinkPEP
        subplot(4,8,28)
%         plot(Tgs, Vgs(:,38), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,38), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,38), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{PEP}')
        xlim([-100 2000])

        % v_sinkPYR
        subplot(4,8,29)
%         plot(Tgs, Vgs(:,39), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,39), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,39), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{PYR}')
        xlim([-100 2000])

        % v_sinkACE
        subplot(4,8,30)
%         plot(Tgs, Vgs(:,40), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,40), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,40), 'color', colorSet(j,:))
            end
            hold on
        end
        title('Sink_{ACE}')
        xlim([-100 2000])

        %recently added (v_mitoNAD)
        subplot(4,8,31)
%         plot(Tgs, Vgs(:,26), 'k')
        for j = 1:npSets
            if j == setup.refParams
                plot(Tgs{j}, Vgs{j}(:,26), 'k','LineWidth',1.2)
            else
                plot(Tgs{j}, Vgs{j}(:,26), 'color', colorSet(j,:))
            end
            hold on
        end
        title('v_{mitoNAD}')
        xlim([-100 2000])
        % set(gcf,'color','w'); suptitle('fig4. Dynamic glucose pulse flux profile. All')
    end
    end

end

% plots being developed
if(setup.plotResultsMode == 50)
    
    if setup.runGSvanHeerden == 1
        % GP insights section 1
        % NAD/NADH balance
        for i = 1
            if setup.runGSvanHeerden == 1
                figure('name','Metabolite pools. Fluxes related. Dynamic glucose pulse profile','units','normalized','outerposition',[0 0 1 1])

                % NADX
                    % NAD
                    subplot(4,4,[1,2,5,6])
        %             plot(Tgs{j},Ygs{j}(:,7),'-','color','black')
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Ygs{j}(:,7),'-','color','black') 
                        else
                            plot(Tgs{j},Ygs{j}(:,7),'color', colorSet(j,:)),
                        end
                        hold on
                    end
                    xlim([-100 340]), legend('NAD')
                    % NADH
                    subplot(4,4,[9,10,13,14])
        %             plot(Tgs{j},Ygs{j}(:,8),'-','color','black')
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Ygs{j}(:,8),'-','color','black') 
                        else
                            plot(Tgs{j},Ygs{j}(:,8),'color', colorSet(j,:)),
                        end
                        hold on
                    end
                    xlim([-100 340]), legend('NADH')

                    % all plots (v_GAPDH, vADH, v_G3PDH, vmitoNADH)
                    subplot(4,4,[3,4,7,8])
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Vgs{j}(:,6),'-','color','magenta','LineWidth',2) 
                            hold on
                            plot(Tgs{j},Vgs{j}(:,8),'-','color','cyan','LineWidth',2) 
                            hold on
                            plot(Tgs{j},Vgs{j}(:,41),'-','color','magenta','LineStyle','--','LineWidth',2) 
                            hold on
                            plot(Tgs{j},Vgs{j}(:,26),'-','color','magenta','LineStyle',':','LineWidth',2) 
                            hold on
                            plot(Tgs{j},+ Vgs{j}(:,6) - Vgs{j}(:,8) + Vgs{j}(:,41) + Vgs{j}(:,26),'-','linewidth',2,'linestyle',':')
                            hold on
                        else
                            plot(Tgs{j},Vgs{j}(:,6))%,'color', colorSet(j,:)), 
                            hold on
                            plot(Tgs{j},Vgs{j}(:,8))%,'color', colorSet(j,:)), 
                            hold on
                            plot(Tgs{j},Vgs{j}(:,41))%,'color', colorSet(j,:)), 
                            hold on
                            plot(Tgs{j},Vgs{j}(:,26))%,'color', colorSet(j,:)), 
                            hold on
                        end
                        hold on
                    end
                    legend('v_{G3PDH}','v_{GAPDH}','v_{ADH}','v_{mitoNADH}')
                    xlim([-100 340])

                    % v_G3PDH
                    subplot(4,4,15)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Vgs{j}(:,6),'-','color','black') 

                        else
                            plot(Tgs{j},Vgs{j}(:,6),'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340]), legend('vG3PDH')
                    % v_GAPDH
                    subplot(4,4,11)
        %             plot(Tgs{j},Vgs{j}(:,8),'-','color','black')
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Vgs{j}(:,8),'-','color','black') 

                        else
                            plot(Tgs{j},Vgs{j}(:,8),'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340]), legend('vGAPDH')
                    % v_ADH
                    subplot(4,4,12)
        %             plot(Tgs{j},Vgs{j}(:,41),'-','color','black')
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Vgs{j}(:,41),'-','color','black') 

                        else
                            plot(Tgs{j},Vgs{j}(:,41),'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340]), legend('vADH')
                    % v_mitoNADH
                    subplot(4,4,16)
        %             plot(Tgs{j},Vgs{j}(:,26),'-','color','black')
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j},Vgs{j}(:,26),'-','color','black') 

                        else
                            plot(Tgs{j},Vgs{j}(:,26),'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340]), legend('vmitoNADH')
                % set(gcf,'color','w'); suptitle('fig44. Metabolite pools. Fluxes related. Dynamic glucose pulse profile')
            end
            % suptitle content
            suptitleContent = {'During glucose pulse, most of the NADH produced by v_{GAPDH} is recycled by fermentation in v_{ADH} (~80%).';...
                'v_{mitoNADH} recycles 18% and only the remaining 2% is taken up by the glycerol branch.'};
            suptitle(suptitleContent);
        end
        % Pi balances
        for i = 1
            if setup.runGSvanHeerden == 1
                figure('units','normalized','outerposition',[0 0 1 1])

                % subplot(6,4,[9,13]) %Pi
                subplot(6,4,[9,13])
                for j = 1:npSets
                    if j == setup.refParams
                        plot(Tgs{j}, Ygs{j}(:,27), 'k','LineWidth',1.2), 
                    else
                        plot(Tgs{j}, Ygs{j}(:,27), 'color', colorSet(j,:)), 
                    end
                    hold on
                end
                xlim([-100 340])
                title('Pi')
                    % subplot(6,4,2) %glycolysis (-gapdh)
                    subplot(6,4,2)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j}, -Vgs{j}(:,8), 'k','LineWidth',1.2), 
                        else
                            plot(Tgs{j}, -Vgs{j}(:,8), 'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340])
                    title('glycolysis (-gapdh)')
                        %subplot(6,4,3) empty
                        subplot(6,4,3)
                        %subplot(6,4,4) %-gapdh
                        subplot(6,4,4)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,8), 'k','color','red'), 
                            else
                                plot(Tgs{j}, Vgs{j}(:,8), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('gapdh')
                    % subplot(6,4,6) %trehalose (+tps1+tps2)
                    subplot(6,4,6)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j}, Vgs{j}(:,21)+Vgs{j}(:,19), 'k','color','black','LineWidth',1.2), 
                        else
                            plot(Tgs{j}, Vgs{j}(:,21)+Vgs{j}(:,19), 'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340])
                    title('trehalose cycle')
                        %subplot(6,4,7) %+tps1+tps2
                        subplot(6,4,7)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,21), 'k','color','blue'),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,19), 'k','color','blue','linestyle',':'),
                                hold on
                            else
                                plot(Tgs{j}, Vgs{j}(:,21), 'color', colorSet(j,:)),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,19), 'color', colorSet(j,:),'linestyle',':'),
                                hold on 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('tps1','tps2')                    
                        %subplot(6,4,8) empty
                        subplot(6,4,8)
                    % subplot(6,4,10) %axp cycle (+atpase-mito)
                    subplot(6,4,10)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j}, Vgs{j}(:,27)-Vgs{j}(:,28), 'k','color','black','LineWidth',1.2), 
                        else
                            plot(Tgs{j}, Vgs{j}(:,27)-Vgs{j}(:,28), 'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340])
                    title('axp cycle')
                        %subplot(6,4,11) %+atpase
                        subplot(6,4,11)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,27), 'k','color','blue'), 
                            else
                                plot(Tgs{j}, Vgs{j}(:,27), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('atpase')
                        %subplot(6,4,12) %-mito
                        subplot(6,4,12)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,28), 'k','color','red'), 
                            else
                                plot(Tgs{j}, Vgs{j}(:,28), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('mito')
                    % subplot(6,4,14) %ixp cycle (+hpt1-isn1)
                    subplot(6,4,14)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j}, Vgs{j}(:,33)-Vgs{j}(:,31), 'k','color','black','LineWidth',1.2), 
                        else
                            plot(Tgs{j}, Vgs{j}(:,33)-Vgs{j}(:,31), 'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340])
                    title('ixp cycle')
                        %subplot(6,4,15) %+hpt1
                        subplot(6,4,15)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,33), 'k','color','blue'), 
                            else
                                plot(Tgs{j}, Vgs{j}(:,33), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('hpt1')
                        %subplot(6,4,16) %-isn1
                        subplot(6,4,16)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,31), 'k','color','red'), 
                            else
                                plot(Tgs{j}, Vgs{j}(:,31), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('isn1')
                    % subplot(6,4,18) %vacuole (+vac import)
                    subplot(6,4,18)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j}, Vgs{j}(:,42), 'k','color','black','LineWidth',1.2), 
                        else
                            plot(Tgs{j}, Vgs{j}(:,42), 'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340])
                    title('vac import')
                        %subplot(6,4,19) %+vac_import
                        subplot(6,4,19)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,42), 'k','color','blue'), 
                            else
                                plot(Tgs{j}, Vgs{j}(:,42), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('vac import')
                        %subplot(6,4,20) empty
                        subplot(6,4,20)
                    % subplot(6,4,22) %sinks (+sinkg6p+sinkp3g+sinkpep-sinkf6p-sinkgap)
                    subplot(6,4,22)
                    for j = 1:npSets
                        if j == setup.refParams
                            plot(Tgs{j}, Vgs{j}(:,34) + Vgs{j}(:,37) + Vgs{j}(:,38) - Vgs{j}(:,35) - Vgs{j}(:,36), 'k','color','black','LineWidth',1.2), 
                        else
                            plot(Tgs{j}, Vgs{j}(:,34) + Vgs{j}(:,37) + Vgs{j}(:,38) - Vgs{j}(:,35) - Vgs{j}(:,36), 'color', colorSet(j,:)), 
                        end
                        hold on
                    end
                    xlim([-100 340])
                    title('sinks overall')
                        % subplot(6,4,23) %sinks (+sinkg6p+sinkp3g+sinkpep)
                        subplot(6,4,23)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,34), 'k','color','blue'),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,37), 'k','color','blue','LineStyle',':'),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,38), 'k','color','blue','LineStyle','--'),
                            else
                                plot(Tgs{j}, Vgs{j}(:,34), 'color', colorSet(j,:)),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,37), 'color', colorSet(j,:),'LineStyle',':'), 
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,38), 'color', colorSet(j,:),'LineStyle','--'), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('sink_{g6p}','sink_{p3g}','sink_{pep}')
                        % subplot(6,4,24) %sinks (-sinkf6p-sinkgap)
                        subplot(6,4,24)
                        for j = 1:npSets
                            if j == setup.refParams
                                plot(Tgs{j}, Vgs{j}(:,35), 'k','color','red'),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,36), 'k','color','red','LineStyle',':'),
                            else
                                plot(Tgs{j}, Vgs{j}(:,35), 'color', colorSet(j,:)),
                                hold on
                                plot(Tgs{j}, Vgs{j}(:,36), 'color', colorSet(j,:)), 
                            end
                            hold on
                        end
                        xlim([-100 340])
                        legend('sink_{f6p}','sink_{gap}')
                % suptitle content
                suptitleContent = {'\color{red}Contribution of ATPase seems to be now the most importat for Pi balance';...
                    'Trehalose and IXP cycle have also an effect, but minor. Uncler Pi accumulation'};
                suptitle(suptitleContent);
            end
        end
        % AXP/IXP balances
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])

            % row 1. Change AXP->IXP
            subplot(3,6,2) %sumAXP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('sum_{AXP}')
            subplot(3,6,3) %sumIXP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('sum_{IXP}')
            subplot(3,6,5) %vAmd1
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,29), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,29), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{Amd1}')
            subplot(3,6,6) %vAde12Ade13
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,30), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,30), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{Ade12Ade13}')

            % row 2. AXP pool internal
            subplot(3,6,7) % ATP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,9), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,9), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('ATP')
            subplot(3,6,8) % ADP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,15), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,15), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('ADP')
            subplot(3,6,9) % AMP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,16), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,16), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('AMP')
            subplot(3,6,10) % v_glycolysis overall
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, - Vgs{j}(:,2) - Vgs{j}(:,4) + Vgs{j}(:,9) + Vgs{j}(:,12), 'k-'),
                    hold on
                    plot(Tgs{j}, - Vgs{j}(:,21), 'r-'),
                else
                    plot(Tgs{j}, - Vgs{j}(:,2) - Vgs{j}(:,4) + Vgs{j}(:,9) + Vgs{j}(:,12), 'color', colorSet(j,:)),
                    hold on
                    plot(Tgs{j}, - Vgs{j}(:,21), 'color', colorSet(j,:)),
                end
                hold on
            end
            xlim([-100 340])
    %         v(9)=  
    %         - v_GLK - v_PFK 
    %         + v_PGK + v_PYK 
    %         - v_TPS1 
    %         title('v^{glycolysis}_{overall}, -GLK-PFK+PGK+PYK\color{red}-TPS1')
            title('v^{glycolysis overall}_{-GLK-PFK+PGK+PYK\color{red}-TPS1}')
            subplot(3,6,11) % v_ATPase
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,27), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,27), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{ATPase}')
            subplot(3,6,12) % v_mito
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,28), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,28), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{mito}')

            % row 3. IXP pool internal
            subplot(3,6,13) % IMP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,28), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,28), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('IMP')
            subplot(3,6,14) % INO
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,29), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,29), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('INO')
            subplot(3,6,15) % HYP
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,30), 'k-'), 
                else
                    plot(Tgs{j}, Ygs{j}(:,30), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('HYP')
            subplot(3,6,16) % v_Isn1
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,31), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,31), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{Isn1}')
            subplot(3,6,17) % v_Pnp1
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,32), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,32), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{Pnp1}')
            subplot(3,6,18) % v_Hpt1
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,33), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,33), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{Hpt1}')

            suptitleText = {'\color{red}Exchange AXP-IXP, internal AXP fluxes, internal IXP fluxes'; 'comment on the transfer between sumAXP and sumIXP'; 'comment on UG first consuming ATP and LG then producing it.'};
            suptitle(suptitleText);
        end

        % GP insights section 2
        % What happens to the trehalose node?
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])

            % subplot(2,4,1) %v_{hxk},v_{pgi}
            subplot(2,4,1) %v_{hxk},v_{pgi}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,2), 'k-'), 
                    hold on
                    plot(Tgs{j}, Vgs{j}(:,3), 'k--'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,2), 'color', colorSet(j,:)), 
                    hold on
                    plot(Tgs{j}, Vgs{j}(:,3), 'color', colorSet(j,:),'LineStyle','--'), 
                end
                hold on
            end
            xlim([-100 340])
            legend('v_{hxk}','v_{pgi}')
            title('v_{hxk},v_{pgi}')
            % subplot(2,4,5) %v_{pgm1},v_{tps1},v_{sinkG6P}
            subplot(2,4,5) %v_{pgm1},v_{tps1},v_{sinkG6P}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,17), 'k-'), 
                    hold on
                    plot(Tgs{j}, Vgs{j}(:,21), 'k--'), 
                    hold on
                    plot(Tgs{j}, Vgs{j}(:,34), 'k:'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,17), 'color', colorSet(j,:)), 
                    hold on
                    plot(Tgs{j}, Vgs{j}(:,21), 'color', colorSet(j,:),'LineStyle','--'), 
                    hold on
                    plot(Tgs{j}, Vgs{j}(:,34), 'color', colorSet(j,:),'LineStyle',':'), 
                end
                hold on
            end
            xlim([-100 340])
            ylim([0 0.1])
            legend('v_{pgm1}','v_{tps1}','v_{sinkG6P}')
            title('v_{pgm1},v_{tps1},v_{sinkG6P}')
            % subplot(2,4,2) %v_{pgi}/v_{hxk}
            subplot(2,4,2) %v_{pgi}/v_{hxk}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,3)./Vgs{j}(:,2), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,3)./Vgs{j}(:,2), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{pgi}/v_{hxk}')
            % subplot(2,4,6) %v_{pgm1}/v_{hxk}
            subplot(2,4,6) %v_{pgm1}/v_{hxk}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,17)./Vgs{j}(:,2), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,17)./Vgs{j}(:,2), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            ylim([0 0.2])
            title('v_{pgm1}/v_{hxk}')
            % subplot(2,4,7) %v_{tps1}/v_{hxk}
            subplot(2,4,7) %v_{tps1}/v_{hxk}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,21)./Vgs{j}(:,2), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,21)./Vgs{j}(:,2), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{tps1}/v_{hxk}')
            % subplot(2,4,8) %v_{sinkg6p}/v_{hxk}
            subplot(2,4,8) %v_{sinkg6p}/v_{hxk}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,34)./Vgs{j}(:,2), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,34)./Vgs{j}(:,2), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkg6p}/v_{hxk}')

            suptitle('What happens to the trehalose node?')
        end
        % What goes into the sinks?
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])

            % subplot(4,5,1) %v_{glt}
            subplot(4,5,1) %v_{glt}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,1), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,1), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{glt}')
            % subplot(4,5,2) %v_{adh}
            subplot(4,5,2) %v_{adh}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,41), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,41), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{adh}')
            % subplot(4,5,4) %v_{sinkG6P}
            subplot(4,5,4) %v_{sinkG6P}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,34), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,34), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkG6P}')
            % subplot(4,5,5) %v_{sinkF6P}
            subplot(4,5,5) %v_{sinkF6P}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,35), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,35), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkF6P}')
            % subplot(4,5,6) %v_{sinkGAP}
            subplot(4,5,6) %v_{sinkGAP}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,36), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,36), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkGAP}')
            % subplot(4,5,7) %v_{sinkP3G}
            subplot(4,5,7) %v_{sinkP3G}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,37), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,37), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkP3G}')
            % subplot(4,5,8) %v_{sinkPEP}
            subplot(4,5,8) %v_{sinkPEP}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,38), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,38), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkPEP}')
            % subplot(4,5,9) %v_{sinkPYR}
            subplot(4,5,9) %v_{sinkPYR}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,39), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,39), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkPYR}')
            % subplot(4,5,10) %v_{sinkACE}
            subplot(4,5,10) %v_{sinkACE}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,40), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,40), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinkACE}')

            % subplot(4,5,[11 12 16 17]) %v_{adh}/v_{glt}
            subplot(4,5,[11 12 16 17]) %v_{adh}/v_{glt}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,41)./Vgs{j}(:,1), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,41)./Vgs{j}(:,1), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{adh}/v_{glt}')

            % subplot(4,5,[14 15 19 20]) %v_{sinksTotal}/v_{glt}
            subplot(4,5,[14 15 19 20]) %v_{sinksTotal}/v_{glt}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, (Vgs{j}(:,34)+Vgs{j}(:,35)+Vgs{j}(:,36)+Vgs{j}(:,37)+Vgs{j}(:,38)+Vgs{j}(:,39)+Vgs{j}(:,40))./Vgs{j}(:,1), 'k-'), 
                else
                    plot(Tgs{j}, (Vgs{j}(:,34)+Vgs{j}(:,35)+Vgs{j}(:,36)+Vgs{j}(:,37)+Vgs{j}(:,38)+Vgs{j}(:,39)+Vgs{j}(:,40))./Vgs{j}(:,1), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{sinksTotal}/v_{glt}')

            suptitle('What goes into the sinks?')

        end
        % Ratio UG_LG
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])
            % subplot(2,3,1) %v_{ald}
            subplot(2,3,1) %v_{ald}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,5), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,5), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{ald}') 
            % subplot(2,3,2) %v_{gpd}
            subplot(2,3,2) %v_{gpd}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,6), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,6), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{gpd}') 
            % subplot(2,3,3) %v_{gapdh}
            subplot(2,3,3) %v_{gapdh}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,8), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,8), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{gapdh}') 
            % subplot(2,3,5) %v_{gpd}/v_{ald}
            subplot(2,3,5) %v_{gpd}/v_{ald}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,6)./Vgs{j}(:,5), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,6)./Vgs{j}(:,5), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{gpd}/v_{ald}')
            % subplot(2,3,6) %v_{gapdh}/v_{ald}
            subplot(2,3,6) %v_{gapdh}/v_{ald}
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,8)./Vgs{j}(:,5), 'k-'), 
                else
                    plot(Tgs{j}, Vgs{j}(:,8)./Vgs{j}(:,5), 'color', colorSet(j,:)), 
                end
                hold on
            end
            xlim([-100 340])
            title('v_{gapdh}/v_{ald}')
        end
        % Any other plot from CA_suarez?
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])
            title('Any other plot from CA_{suarez}?')
        end
        % Which is the main stopper/regulator?
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])
            title('Which is the main stopper/regulator?')
        end
        % 'How cofactors affect glycolysis? What if they are clamped? (benefit of integrated view)'
        for i = 1
            figure('units','normalized','outerposition',[0 0 1 1])
            suptitleName = {'How cofactors affect glycolysis? What if they are clamped? What is their effect on glycolysis? (benefit of integrated view)';...
                            'clamp tre: where does it affect?';...
                            'clamp pi: where does it affect?';...
                            'clamp nadx: where does it affect?';...
                            'clamp axp: where does it affect?';...
                            'clamp ixp: where does it affect?';...
                            '\color{red}Add SS study here. To show the correlation between glycerol';...
                            'and NADH metabolism and how cofactors (like AXP) also regulate in SS'};
            title(suptitleName);
        end
    end
    
    
    if setup.runSScanelas == 1
        % SS insights. Flux distribution
        % 1.1. SS simulations. Compartments vs main glycolysis
        for i = 1
            mu = dprofile{1};
            flux_GLT = exp_v_GLT{1}';
            flux_vADH = exp_v_ADH{1}'/2;
            flux_vsinkPYR = exp_v_sinkPYR{1}'/2;
            flux_vsinkACE = exp_v_sinkACE{1}'/2;
            flux_vglycbranch = exp_v_G3PDH{1}';
            flux_vsinkG6P = exp_v_sinkG6P{1}';
            flux_vsinkF6P = exp_v_sinkF6P{1}';
            flux_vsinkGAP = exp_v_sinkGAP{1}';
            flux_vsinkP3G = exp_v_sinkP3G{1}'/2;
            flux_vsinkPEP = exp_v_sinkPEP{1}'/2;
            totalEntry = flux_GLT + flux_vsinkF6P + flux_vsinkGAP;
            
            percent_sinkpyr = flux_vsinkPYR./totalEntry;
            percent_vADH = flux_vADH./totalEntry;
            percent_sinkACE = flux_vsinkACE./totalEntry;
            percent_vglycbranch = flux_vglycbranch./totalEntry;
            percent_sinkG6P = flux_vsinkG6P./totalEntry;
            percent_sinkP3G = flux_vsinkP3G./totalEntry;
            percent_sinkPEP = flux_vsinkPEP./totalEntry;

            figure('units','normalized','outerposition',[0 0 1 1])
            
            subplot(4,1,[2,3,4])
            h1 = area(mu,[percent_sinkpyr',percent_vADH',percent_sinkACE',percent_vglycbranch',percent_sinkG6P',percent_sinkP3G',percent_sinkPEP'],'LineStyle','none');
            h1(5).FaceColor = h1(3).FaceColor;
            h1(6).FaceColor = [0.9290+(1-0.9290)*1/3    0.6940+(1-0.6940)*1/3    0.1250+(1-0.1250)*1/3];
            h1(7).FaceColor = [0.9290+(1-0.9290)*2/3    0.6940+(1-0.6940)*2/3    0.1250+(1-0.1250)*2/3];
            h1(3).FaceColor = [0.8500+(1-0.8500)*1/3    0.3250+(1-0.3250)*1/3    0.0980+(1-0.0980)*1/3];
            h1(4).FaceColor = [0.8500+(1-0.8500)*2/3    0.3250+(1-0.3250)*2/3    0.0980+(1-0.0980)*2/3];
            xlim([mu(1) mu(end)])
            ylabel('% of GLT entry flux')
            legend('sinkPYR','vADH','sinkACE','vglycbranch','sinkG6P','sinkP3G','sinkPEP','location','southoutside','orientation','horizontal')
            
            subplot(4,1,1)
            area(mu, flux_GLT,'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
            ylabel('amount [??]')
            xlabel('dilution rate [h^{-1}]')
            xlim([mu(1) mu(end)])
            %?? add colors per groups
            %?? add %trehalose
            %?? then correct part of sinkG6P that comes back in sinkF6P or sinkGAP
            %?? add safechek at height 0
            suptitleName = {'SS simulations. Compartments vs main glycolysis';'GLT flux vs the rest of the system'};
            suptitle(suptitleName);
        end
        % 1.2.1. SS simulations. Cofactors: AXP
        for i = 1
            ATPglycolysis = - exp_v_HK{1} - exp_v_PFK{1} + exp_v_PGK{1} + exp_v_PYK{1}; % calc ATP produced in glycolysis
            ATPmitochondria = + exp_v_mito{1}; % calc ATP produced in mitochondria
            ATPtotal = ATPglycolysis + ATPmitochondria; % calc ATP produced total
            flux_glycolysis = ATPglycolysis./ATPtotal;
            flux_mitochondria = ATPmitochondria./ATPtotal;
            figure('units','normalized','outerposition',[0 0 1 1])
            
            subplot(4,1,[2,3,4])
            area(mu,[flux_glycolysis,flux_mitochondria],'LineStyle','none')
            xlim([mu(1) mu(end)])
            ylabel('% of ATP synthesis')
            legend('glycolysis','mitochondria','location','southoutside','orientation','horizontal')
            
            subplot(4,1,1)
            area(mu, ATPtotal','FaceColor',[0.9 0.9 0.9],'LineStyle','none')
            ylabel('amount [??]')
            xlabel('dilution rate [h^{-1}]')
            xlim([mu(1) mu(end)])
            
            suptitle('SS simulations. Cofactors: ATP synthesized and origin')
        end
        % 1.2.2. SS simulations. Cofactors: NADX
        for i = 1
            rate_GAPDH = exp_v_GAPDH{1};
            rate_G3PDH = exp_v_G3PDH{1};
            rate_ADH = exp_v_ADH{1};
            rate_mitoNADH = exp_v_mitoNADH{1};
            flux_G3PDH = rate_G3PDH ./ rate_GAPDH;
            flux_ADH = rate_ADH ./ rate_GAPDH;
            flux_mitoNADH = rate_mitoNADH ./ rate_GAPDH;
            
            figure('units','normalized','outerposition',[0 0 1 1])
            
            subplot(4,1,[2,3,4])
            h1 = area(mu,[flux_mitoNADH,flux_ADH,flux_G3PDH],'LineStyle','none');
            h1(3).FaceColor = [0.9000 0.6500 0.5000];
            xlim([mu(1) mu(end)])
            ylim([0 1])
            ylabel('% of NADX synthesis')
            legend('flux_{mitoNADH}','flux_{ADH}','flux_{G3PDH}','location','southoutside','orientation','horizontal')
            
            subplot(4,1,1)
            area(mu, rate_GAPDH','FaceColor',[0.9 0.9 0.9],'LineStyle','none')
            ylabel('amount [??]')
            xlabel('dilution rate [h^{-1}]')
            xlim([mu(1) mu(end)])
            
            suptitle('SS simulations. Cofactors: NADX metabolism')
        end
    end
end

%%
% plots physiology (for I = 1)
if(setup.plotResultsMode == 70)

%     figure
%     cell23={{['-g']},{['-g']},{['-g']};{['-g']},{['-g']},{['-g']}};
%     C = cell23;
%     [h,labelfontsize] = subplotplus(C);

    % fig21. Physiology. Steady state profile. All
    for i = 1
    if setup.runSScanelas == 1
% % % %         figure('name','Physiology. Steady state profile. All','units','normalized','outerposition',[0 0 1 1])
        % qS
%         subplot(2,3,1)
        
        subidx = 1;
        set(gcf,'CurrentAxes',h(subidx));
        
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'ko-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',5)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GLT, 'kd',...
            'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',5)
        end
%         title('q_{S}')
        xlim([0 0.4])
%         ylim([0 4])
        ylim([0 2])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
            xt = 0.05;
            yt = 1.75;
            str = 'a';
            text(xt,yt,str,'FontSize',20)
        % qEtoh
%         subplot(2,3,4)
        
        subidx = 4;
        set(gcf,'CurrentAxes',h(subidx));
        
        for j = 1:npSets
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'ko-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',5)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH, 'kd',...
            'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',5)
        end
%         title('q_{Etoh}')
        xlim([0 0.4])
%         ylim([0 4])
        ylim([0 2])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
            xt = 0.05;
            yt = 1.75;
            str = 'd';
            text(xt,yt,str,'FontSize',20)
        % qCO2
%         subplot(2,3,2)
        
        subidx = 5;
        set(gcf,'CurrentAxes',h(subidx));
        
        for j = 1:npSets
            if j == setup.refParams
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, ...
                    'ko-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',5)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'color', colorSet(j,:))
            end
        hold on
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qCO2, 'kd',...
            'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',5)
        end
%         title('q_{CO2}')
        xlim([0 0.4])
        ylim([0 4])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
            xt = 0.05;
            yt = 3.5;
            str = 'e';
            text(xt,yt,str,'FontSize',20)
        % qO2
%         subplot(2,3,5)
        
        subidx = 2;
        set(gcf,'CurrentAxes',h(subidx));
        
        for j = 1:npSets
            if j == setup.refParams
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'k','LineWidth',1.2)
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', ...
                    'ko-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',5)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'k','LineWidth',1.2)
            end
        hold on
        plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qO2,'kd',...
            'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',5)
        end
%         title('q_{O2}')
        xlim([0 0.4])
%         ylim([0 4])
        ylim([0 2])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
            xt = 0.05;
            yt = 1.75;
            str = 'b';
            text(xt,yt,str,'FontSize',20)
% % % %      
        % calc qmito and qATPase
        for k = 27
            qATPase = zeros(1,length(dprofile{1}));
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{1}.ss2can{o}(end,k);
                qmito(o) = Vss{1}.ss2can{o}(end,k+1);
            end
        end
        % qmito
%         subplot(2,3,3)
        
        subidx = 3;
        set(gcf,'CurrentAxes',h(subidx));
        
        for j = 1:npSets
            % calc qmito and qATPase
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qmito(o) = Vss{j}.ss2can{o}(end,k+1);
            end            
            % sims            
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), qmito, 'ko-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',5)
            else
                plot(dprofile{j}(1:q(1)), qmito, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D,(exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]'/1, 'b+')
        plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.sinkPYR * 3 + 1 * canelas_SS.mtlD.v_PDC)*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]/1,...
            'kd','MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',5)
        end
%         title('q_{mito}')
        xlim([0 0.4])
%         ylim([0 4])
        ylim([0 2])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
            xt = 0.05;
            yt = 1.75;
            str = 'c';
            text(xt,yt,str,'FontSize',20)
        
        % qATPase
%         subplot(2,3,6)
        
        subidx = 6;
        set(gcf,'CurrentAxes',h(subidx));
        
        for j = 1:npSets
% % % %                
            qATPase = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{j}.ss2can{o}(end,k);
            end            
% % % %             
            if j == setup.refParams
                plot(dprofile{j}(1:q(1)), qATPase, 'ko-','LineWidth',1,'MarkerFaceColor','white','MarkerSize',5)
            else
                plot(dprofile{j}(1:q(1)), qATPase, 'color', colorSet(j,:))
            end
            hold on
%             plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/1.7/3.6, 'b+')
            plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/2.0/3.6,...
                 'kd','MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',5)
% % % %             hold on
% % % %             plot(canelas_SS.mtlD.D,(40*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'g+')
        end
%         title('q_{ATPase}')
        xlim([0 0.4])
        ylim([0 4])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
            xt = 0.05;
            yt = 3.5;
            str = 'f';
            text(xt,yt,str,'FontSize',20)
        
        % set(gcf,'color','w'); suptitle('fig21. Physiology. Steady state profile. All')
%         set(gcf,'color',backColor);
        set(gcf,'color','w');
    end
    end
    
    %
    for i = 1
        if setup.runGSvanHeerden == 1
            
            figure
            cell32 = {{['-g']},{['-g']};{['-g']},{['-g']};{['-g']},{['-g']}};
            C = cell32;
            [h,labelfontsize] = subplotplus(C);
            
            % ATP balance concentration
            subidx = 1;
            set(gcf,'CurrentAxes',h(subidx));
            plot(Tgs{j}, Vgs{j}(:,21), 'k','LineWidth',1.2)
            xlim([-100 340])
            legend('vTPS1 concentration')
            
            % ATP balance concentration
            subidx = 2;
            set(gcf,'CurrentAxes',h(subidx));
            yyaxis right
            plot(Tgs{j}, Ygs{j}(:,9), 'k','LineWidth',1.2)
            xlim([-100 340])
            grid off
            legend('ATP concentration')
            yyaxis left
            set(gca,'YTickLabel',[]);
            
            % Pi balance concentration
            subidx = 4;
            set(gcf,'CurrentAxes',h(subidx));
            yyaxis right
            plot(Tgs{j}, Ygs{j}(:,27), 'k','LineWidth',1.2)
            xlim([-100 340])
            ylim([0 30])
            grid off
            yyaxis left
            set(gca,'YTickLabel',[]);
            legend('Pi concentration')

            % Pi balance contribution of rates
            subidx = 3;
            set(gcf,'CurrentAxes',h(subidx));
            
            plot(Tgs{j}, -Vgs{j}(:,8), 'k','LineWidth',1.2), 
            hold on
            plot(Tgs{j}, Vgs{j}(:,27)-Vgs{j}(:,28)+Vgs{j}(:,33)-Vgs{j}(:,31), 'k:','LineWidth',1.2), 
            hold on
% % % %             plot(Tgs{j}, Vgs{j}(:,33)-Vgs{j}(:,31), 'k','color','red','LineWidth',1.2),
% % % %             hold on
            plot(Tgs{j}, Vgs{j}(:,42), 'k--','LineWidth',1.2),
            
            xlim([-100 340])
            ylim([-0.5 0.5])
            legend('glycolysis uptake (gapdh)','atp hydrolysis, synthesis and salvage','vac import')
            
            % AXP/IXP balances
            subidx = 6;
            set(gcf,'CurrentAxes',h(subidx));
%             set(gca,'YTickLabel',[],'defaultAxesColorOrder',[[1 1 1]; [1 1 1]]);  
            yyaxis right
            plot(Tgs{j}, Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16), 'k-','LineWidth',1.2), 
            hold on
            plot(Tgs{j}, Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30), 'k--','LineWidth',1.2), 
            xlim([-100 340])
            grid off
            legend('sum_{AXP}','sum_{IXP}')
            yyaxis left
            set(gca,'YTickLabel',[]);
%             set(gca,'YTickLabel',[],'defaultAxesColorOrder',[[1 1 1]; [1 1 1]]);           
            
            % NAD/NADH balance
            subidx = 5;
            set(gcf,'CurrentAxes',h(subidx));
            for j = 1:npSets
                plot(Tgs{j},Vgs{j}(:,8),'-','color','black','LineWidth',1.2) 
                hold on
                plot(Tgs{j},Vgs{j}(:,41),'-','color','black','LineStyle','--','LineWidth',2) 
                hold on
                plot(Tgs{j},Vgs{j}(:,26),'-','color','black','LineStyle',':','LineWidth',2) 
            end
            legend('NADH production','NADH consumption (fermentation)','NADH consumption (mitochondria)')
            xlim([-100 340])
            
            
            % %$ %% % % % 
            figure
% % % %             cell32 = {{['-g']},{['-g']};{['-g']},{['-g']};{['-g']},{['-g']}};
% % % %             C = cell32;
            cell13 = {{['-g']},{['-g']},{['-g']}};
            C = cell13;
            [h,labelfontsize] = subplotplus(C);
            
            % ATP balance concentration
            subidx = 1;
            set(gcf,'CurrentAxes',h(subidx));
            plot(Tgs{j}, Vgs{j}(:,21), 'k','LineWidth',1.2)
            hold on
            plot(setup.simResults2{1}.gs.T, setup.simResults2{1}.gs.V(:,21), 'color','red','LineWidth',2)
            xlim([-100 340])
            legend('vTPS1 concentration')
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);
                xt = -90;
                yt = 0.023;
                str = 'a';
                text(xt,yt,str,'FontSize',20)
            
            % ATP balance concentration
            subidx = 2;
            set(gcf,'CurrentAxes',h(subidx));
            yyaxis right
            plot(Tgs{j}, Ygs{j}(:,9), 'k','LineWidth',1.2)
            hold on
            plot(setup.simResults2{1}.gs.T, setup.simResults2{1}.gs.Y(:,9), 'r-','LineWidth',1.2)
            xlim([-100 340])
            grid off
            legend('ATP concentration')
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);
                xt = -90;
                yt = 3.05;
                str = 'b';
%                 str = '\color{red}b missing the mutant strain';
                text(xt,yt,str,'FontSize',20)
            yyaxis left
            set(gca,'YTickLabel',[]);
            
% % % %             % Pi balance concentration
% % % %             subidx = 4;
% % % %             set(gcf,'CurrentAxes',h(subidx));
% % % %             yyaxis right
% % % %             plot(Tgs{j}, Ygs{j}(:,27), 'k','LineWidth',1.2)
% % % %             xlim([-100 340])
% % % %             ylim([0 30])
% % % %                 set(gca,'YTickLabel',[]);
% % % %                 set(gca,'XTickLabel',[]);
% % % %                 xt = -90;
% % % %                 yt = 27.5;
% % % %                 str = 'd';
% % % %                 text(xt,yt,str,'FontSize',20)
% % % %             grid off
% % % %             yyaxis left
% % % %             set(gca,'YTickLabel',[]);
% % % %             legend('Pi concentration')
% % % % 
% % % %             % Pi balance contribution of rates
% % % %             subidx = 3;
% % % %             set(gcf,'CurrentAxes',h(subidx));
% % % %             
% % % %             plot(Tgs{j}, -Vgs{j}(:,8), 'k','LineWidth',1.2), 
% % % %             hold on
% % % %             plot(Tgs{j}, Vgs{j}(:,27)-Vgs{j}(:,28)+Vgs{j}(:,33)-Vgs{j}(:,31), 'k:','LineWidth',1.2), 
% % % %             hold on
% % % % % % % %             plot(Tgs{j}, Vgs{j}(:,33)-Vgs{j}(:,31), 'k','color','red','LineWidth',1.2),
% % % % % % % %             hold on
% % % % %             plot(Tgs{j}, Vgs{j}(:,42), 'k--','LineWidth',1.2),
% % % %             Tsel = Tgs{j}([1:3001,3011:end]);
% % % %             Vsel = Vgs{j}([1:3001,3011:end],42);
% % % %             plot(Tsel, Vsel, 'k--','LineWidth',1.2),
% % % %                 set(gca,'YTickLabel',[]);
% % % %                 set(gca,'XTickLabel',[]);
% % % %                 xt = -90;
% % % %                 yt = 0.4;
% % % %                 str = 'c';
% % % %                 text(xt,yt,str,'FontSize',20)
% % % %             
% % % %             xlim([-100 340])
% % % %             ylim([-0.5 0.5])
% % % %             legend('glycolysis uptake (gapdh)','atp hydrolysis, synthesis and salvage','vac import')
            
            % AXP/IXP balances
            subidx = 3;
            set(gcf,'CurrentAxes',h(subidx));
%             set(gca,'YTickLabel',[],'defaultAxesColorOrder',[[1 1 1]; [1 1 1]]);  
            yyaxis right
            plot(Tgs{j}, Ygs{j}(:,9) + Ygs{j}(:,15) + Ygs{j}(:,16), 'k-','LineWidth',1.2), 
            hold on
            plot(data.time_nucleotides, data.AXP + (4.5-2.7), 'o','MarkerSize',7.5,'MarkerEdgeColor','k','MarkerFaceColor','none'),
            hold on
            plot(Tgs{j}, Ygs{j}(:,28) + Ygs{j}(:,29) + Ygs{j}(:,30), 'k--','LineWidth',1.2), 
            xlim([-100 340])
            grid off
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);
                xt = -90;
                yt = 4.25;
                str = 'c';
                text(xt,yt,str,'FontSize',20)
            legend('sum_{AXP,sim}','sum_{AXP,exp}','sum_{IXP}')
            yyaxis left
            set(gca,'YTickLabel',[]);
%             set(gca,'YTickLabel',[],'defaultAxesColorOrder',[[1 1 1]; [1 1 1]]);           
            
% % % %             % NAD/NADH balance
% % % %             subidx = 5;
% % % %             set(gcf,'CurrentAxes',h(subidx));
% % % %             for j = 1:npSets
% % % %                 plot(Tgs{j},Vgs{j}(:,8),'-','color','black','LineWidth',1.2) 
% % % %                 hold on
% % % %                 plot(Tgs{j},Vgs{j}(:,41),'-','color','black','LineStyle','--','LineWidth',2) 
% % % %                 hold on
% % % %                 plot(Tgs{j},Vgs{j}(:,26),'-','color','black','LineStyle',':','LineWidth',2) 
% % % %             end
% % % %             legend('NADH production','NADH consumption (fermentation)','NADH consumption (mitochondria)')
% % % %             xlim([-100 340])
% % % %                 set(gca,'YTickLabel',[]);
% % % %                 set(gca,'XTickLabel',[]);
% % % %                 xt = -90;
% % % %                 yt = 0.45;
% % % %                 str = 'e';
% % % %                 text(xt,yt,str,'FontSize',20)
                
        set(gcf,'color','w');
        
        end
    end    
    
end

%%
% plots profiles
if(setup.plotResultsMode == 80)
    %fig3. Dynamic glucose pulse metabolite profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
%         figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        figure(8001)
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.metabolites.onoff(k) == 1
                plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            xlim(setup.vHplotRange)
% % % %             xlim([-3000 -2500])
% % % %             xlim([-100 340])
% % % %             xlim([-100 3400])
            if k == 25
                xlabel('time [s]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig3. Dynamic glucose pulse metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
    %fig4. Dynamic glucose pulse flux profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
%         figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        figure(8002)
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,k), 'k','LineWidth',1.2)
                else
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.fluxes.onoff(k) == 1
                plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            xlim(setup.vHplotRange)
% % % %             xlim([-3000 -2500])
% % % %             xlim([-100 340])
% % % %             xlim([-100 3400])
            if k == 37
                xlabel('time [s]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end       
        suptitle('fig4. Dynamic glucose pulse flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
end

%% num 81
% plots profiles
if(setup.plotResultsMode == 81)

    % fig1. Steady state metabolite profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Yss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                elseif j == 2
                    plot(dprofile{j}, tempVal, 'color', [0.75 0.75 0.75], 'LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.metabolites.onoff(k) == 1
                errorbar(canelas_SS.ordered.dprof{k},canelas_SS.ordered.metabolites.data{k},canelas_SS.ordered.metabolites.std{k},'r.')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            if k == 25
                xlabel('dilution rate [h^{-1}]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig1. Steady state metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    % fig2. Steady state flux profile. All
    for i = 1
    if setup.runSScanelas == 1
        figure('name','Steady state flux profile. All','units','normalized','outerposition',[0 0 1 1])
    
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Vss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                elseif j == 2
                    plot(dprofile{j}, tempVal, 'color', [0.75 0.75 0.75], 'LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.fluxes.onoff(k) == 1
                plot(canelas_SS.ordered.dprof{end},canelas_SS.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            if k == 37
                xlabel('dilution rate [h^{-1}]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end
       
        suptitle('fig2. Steady state flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end

    %fig3. Dynamic glucose pulse metabolite profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        for k = 1:length(legenda.metabolites)
            subplot(4,8,k)
            % plot simulations
            for j = 1:npSets
%                 if j == setup.refParams
                if j == 1
                    plot(Tgs{j}, Ygs{j}(:,k), 'k','LineWidth',1.2)
                elseif j == 2
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', [0.75 0.75 0.75], 'LineWidth',1.2)
                else
%                     plot(Tgs{j}, Ygs{j}(:,k), 'color', 'blue')
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
                end
%                 disp(j);
                hold on
            end
            % plot experimental data
            if data.ordered.metabolites.onoff(k) == 1
                plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
            end
            titleName = erase(legenda.metabolites{k},", [mM]");
            title(titleName)%fontsize{16}
            xlim([-100 340])
% % % %             xlim([-3000 340])
% % % %             xlim([-100 3400])
% % % %             xlim([-3000 -2500])
            
            if k == 25
                xlabel('time [s]')
                ylabel('concentration [mM]')
            end
        end
        suptitle('fig3. Dynamic glucose pulse metabolite profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
    
    %fig4. Dynamic glucose pulse flux profile. All
    for i = 1
    if setup.runGSvanHeerden == 1
        figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
        
        for k = 1:length(legenda.fluxes)
            subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
%                 if j == setup.refParams
                if j == 1
                    plot(Tgs{j}, Vgs{j}(:,k), 'k','LineWidth',1.2)
                elseif j == 2
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', [0.75 0.75 0.75], 'LineWidth',1.2)
                else
%                     plot(Tgs{j}, Vgs{j}(:,k), 'color', 'blue')
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.fluxes.onoff(k) == 1
                plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},'r+')
            end
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
            title(titleName)%fontsize{16}
            xlim([-100 340])
% % % %             xlim([-3000 340])
% % % %             xlim([-100 3400])
% % % %             xlim([-3000 -2500])
            
            if k == 37
                xlabel('time [s]')
                ylabel('reaction rate [mM s^{-1}]')
            end
        end  
        
        legendaNames = cell(1,length(npSets));
        for j = 1:npSets
            legendaNames{j} = num2str(j);
        end
        legend(legendaNames)  
        
        suptitle('fig4. Dynamic glucose pulse flux profile. All')
        % set(gcf,'color','w');
        set(gcf,'color',backColor);
    end
    end
    
end


% cleaning up things
numFigs = setup.plotResultsMode; % add here the resulting validation I would say
end

