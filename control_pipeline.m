% % Control pipeline
%  The aim of this file is to check rerunning codes from start works well
%  after simplification of the pipeline.
check_confirm = 1;
% Please, run this code in sections, as figure numbering is important for
% some plots


%% f1
clearvars -except check_confirm
close all
F1_reference_simulations
if check_confirm == 1
    fh_current = figure(1);
    fh_current_gp = figure(3);
    openfig('F1_SSmets.fig');
    fh_safecopy = figure(6);
    openfig('F1_GPmets.fig');
    fh_safecopy_gp = figure(7);
    if((sum(fh_current.Children(end).Children.YData == fh_safecopy.Children(end).Children.YData) == 8)&&(sum(fh_safecopy_gp.Children(end).Children.YData) == sum(fh_current_gp.Children(end).Children.YData)))
        disp('F1 reference simulations are reproduced properly.')
    else
        disp('F1 reference simulations are NOT reproduced properly.')
    end
    close(6)
    close(7)
end


%% f2
clearvars -except check_confirm
close all
F2_physio
if check_confirm == 1
    fh_current = figure(10002);
    openfig('F2_energetics.fig');
    fh_safecopy = figure(7);
    if(sum(fh_current.Children(2).Children(end-1).YData == fh_safecopy.Children(2).Children(end-1).YData) == 8)
        disp('F2 physiology is reproduced properly.')
    else
        disp('F2 physiology is NOT reproduced properly.')
    end
    close(7)
end


%% f3
clearvars -except check_confirm
close all
F3_par_est
if check_confirm == 1
    fh_current = figure(7);
    openfig('F3_fit_GP_PEP.fig');
    fh_safecopy = figure(8);
    if(sum(fh_current.Children(2).Children(4).YData) == sum(fh_safecopy.Children(2).Children(4).YData))
        disp('F3 parameter estimation is reproduced properly.')
    else
        disp('F3 parameter estimation is NOT reproduced properly.')
    end
    close(8)
end


%% f4
clearvars -except check_confirm
close all
F4_data
if check_confirm == 1
    fh_current = figure(6002);
    openfig('F4_estimation.fig');
    fh_safecopy = figure(7);
    if(sum(abs(fh_current.Children(2).Children(end-1).XData)) == sum(abs(fh_safecopy.Children(2).Children(end-1).XData)))
        disp('F4 data is reproduced properly.')
    else
        disp('F4 data is NOT reproduced properly.')
    end
    close(7)
end


%% f5
clearvars -except check_confirm
close all
F5_robustness
if check_confirm == 1
    fh_current = figure(2001);
    openfig('F5_robust_ss_mets.fig');
    fh_safecopy = figure(7);
%     temp_id = round(32 * rand);
    if(fh_current.Children(22).Children(5).CData(60,72) == fh_safecopy.Children(22).Children(5).CData(60,72))
        disp('F5 robustness is reproduced properly.')
    else
        disp('F5 robustness is NOT reproduced properly.')
    end
    close(7)
end


%% f6
clearvars -except check_confirm
close all
F6_lit_param
if check_confirm == 1
    fh_current = figure(100);
    openfig('F6_lit_param_sims.fig');
    fh_safecopy = figure(1);
    if(sum(abs(fh_current.Children(5).Children(2).YData)) == sum(abs(fh_safecopy.Children(5).Children(2).YData)))
        disp('F6 literature parameter simulations are reproduced properly.')
    else
        disp('F6 literature parameter simulations are NOT reproduced properly.')
    end
    close(1)
end


%% f7
clearvars -except check_confirm
close all
F7_hxt_deppend
if check_confirm == 1
    fh_current = figure(5001);
    openfig('F7_glt_mu_dependency.fig');
    fh_safecopy = figure(7);
    if(sum(abs(fh_current.Children(2).Children.YData)) == sum(abs(fh_safecopy.Children(2).Children.YData)))
        disp('F7 HXT dilution rate dependency is reproduced properly.')
    else
        disp('F7 HXT dilution rate dependency is NOT reproduced properly.')
    end
    close(7)
end


%% f8
clearvars -except check_confirm
close all
F8_modules
if check_confirm == 1
    fh_current = figure(7001);
    openfig('F8_sinks.fig');
    fh_safecopy = figure(3);
    if(sum(fh_current.Children(6).Children(2).YData == fh_safecopy.Children(6).Children(2).YData) == 8)
        disp('F8 module simulations are reproduced properly.')
    else
        disp('F8 module simulations are  NOT reproduced properly.')
    end
    close(3)
end





