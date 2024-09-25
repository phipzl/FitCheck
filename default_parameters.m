% default_parameters.m

%% General Settings
reload_data = 1;

%% Paths
lab_dir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab';
sim_dir = 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/Ala_GSH_Simulations';
process_results_dir = fullfile(lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/Results_v2/');
output_dir = fullfile(lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/FittingSimulationResults');
batch_dir = fullfile(output_dir, 'BatchDir');
lcmodel_program_path = '/usr/local/lcmodel/bin/lcmodel';
control_info_file = fullfile(sim_dir, 'ControlFiles/LCModel_Control_GH_2020_Pat_sMM_ForRevision.m');
basis_sets_dir = fullfile(lab_dir, 'Basis_Sets/GH_FID_Basis_2019/jmrui');
basis_set = fullfile(lab_dir, 'Basis_Sets/GH_FID_Basis_2019/Version_1_1_1/fid_0.000000ms.basis');

%% Metabolite Names and Basis Files
mets = {
    'NAA', 'NAAG', ...
    'Cr', 'PCr', ...
    'GPC', 'PCh', ...
    'Glu', 'Gln', ...
    'Gly', 'Ins', ...
    'GSH', 'Ala', ...
    'Tau', 'Ser', 'Cys'
};

met_files = {
    'NAA_damp55.txt', 'NAAG.txt', 'Cr_Damp190.txt', 'PCr_Damp190.txt', 'GPC.txt', 'PCh.txt', ...
    'glu.txt', 'Gln_Damp130.txt', 'glycine.txt', 'Ins.txt', 'GSH_Damp100.txt', 'Ala_fixed.txt', ...
    'Tau.txt', 'Ser.txt', 'Cys.txt'
};

%% Metabolic Composition
% Initialize default levels
default_values = [... 
    90,  30, ... % NAA+NAAG
    40,  20, ... % Cr+PCr
    18,   7, ... % GPC+PCh
    25,  25, ... % Glu, Gln
     5, 100, ... % Gly, Ins
    30,  30, ... % GSH, Ala
    0, 0, 0]; % Tau, Ser, Cys

% Initialize met_lvls_default using 'mets' and 'default_values'
met_lvls_default = cell2struct(num2cell(default_values), mets, 2);

% Variable Levels
mets_variable = {'Ala', 'GSH'};  % Metabolites to vary
met_lvls_variable = {
    [10, 30, 50];  % Levels for 'Ala'
    [10, 30, 50];  % Levels for 'GSH'
};

%% Tissue Types and Parameters
% noise levels. 
% 50: low noise, 55: okay, relatively low noise, 58: still okay
% 59: SNR 14, 60: complete mess
% filter levels
% 1: does almost nothing, % 2: looks okay
% 2.5: wrong fits with low CRLBs, but spectrum looks ok
% 3: completely broken
noise_dbs = 0; % 59:0.25:60; %4e3:4e3:12.2e4;
filter_levels = 1; %1.5:0.5:2.5;


%% Simulation Parameters
dursim = 1.092157;     % Simulation duration
durmeas = 0.3456;      % Measurement duration

%% MetaInfo Parameters
MetaInfo.DimNames = {'X'};
MetaInfo.Dimt1 = 4;
MetaInfo.pat_name = 'My_Simulation';
MetaInfo.LarmorFreq = 297.2232 * 1E6;  % [Hz]
dwt = 3.2E-4;                          % Dwell time in seconds
MetaInfo.dwelltime = dwt * 1E9;        % [ns] for LCModel

%% CPU Cores
CPU_cores = 1;

