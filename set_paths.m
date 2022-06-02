% % set_paths.m
% Setting the working path.
folder=cd;
temp_path = [folder,'\experimental_data'];
    addpath(genpath(temp_path)); clear temp_path
temp_path = [folder,'\model'];
    addpath(genpath(temp_path)); clear temp_path
temp_path = [folder,'\parameter_estimation'];
    addpath(genpath(temp_path)); clear temp_path
temp_path = [folder,'\manuscript_results'];
    addpath(genpath(temp_path)); clear temp_path
temp_path = [folder,'\simulations'];
    addpath(genpath(temp_path)); clear temp_path
temp_path = [folder,'\visualization'];
    addpath(genpath(temp_path)); clear temp_path

