% default_parameters.m

%% General Settings
MetaInfo.reload_data = 0; % Set this to 1 to always reload the spectral data

%% Paths
Paths = struct();
Paths.lab_dir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab';
Paths.sim_dir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/LCM_FitCheck';
Paths.process_results_dir = fullfile(Paths.lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/Results_v2/');
Paths.out_dir = fullfile(Paths.lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/FittingSimulationResults');
Paths.batchdir = fullfile(Paths.out_dir, 'BatchDir');
Paths.LCM_ProgramPath = '/usr/local/lcmodel/bin/lcmodel';
Paths.ControlInfo = fullfile(Paths.sim_dir, 'ControlFiles/LCModel_Control_GH_2020_Pat_sMM_ForRevision.m');
Paths.basis_sets_dir = fullfile(Paths.sim_dir, 'BasisSet1.3ms');
Paths.basis_file = {fullfile(Paths.lab_dir, 'Basis_Sets/GH_FID_Basis_2019/Version_1_1_1/fid_0.000000ms.basis')};

%% Write_LCM_files parameters (incl. MetaInfo)

MetaInfo.DimNames = {'X', 'Y', 'Z'};
MetaInfo.Dimt1 = 4;
MetaInfo.pat_name = 'My_Simulation';
MetaInfo.LarmorFreq = 297.2232 * 1E6;  % [Hz]
MetaInfo.dwelltime_seconds = 3.4E-4;                          % Dwell time in seconds
MetaInfo.dwelltime = MetaInfo.dwelltime_seconds * 1E9;        % [ns] for LCModel

%% Simulation Parameters
MetaInfo.SimDuration  = 1.092157;     % dwt * (number of FID points - 1)
MetaInfo.MeasDuration = 0.345600;     % ADC duration

%% Metabolite Names and Basis Files
Metabolite.names = {
    'NAA', 'NAAG', ...
    'Cr', 'PCr', ...
    'GPC', 'PCh', ...
    'Glu', 'Gln', ...
    'Gly', 'Ins', ...
    'GSH', 'Ala', ...
    'Tau', 'Ser', 'Cys'
};

Metabolite.files = {
    'NAA_damp55.txt', 'NAAG.txt', 'Cr_Damp190.txt', 'PCr_Damp190.txt', 'GPC.txt', 'PCh.txt', ...
    'glu.txt', 'Gln_Damp130.txt', 'glycine.txt', 'Ins.txt', 'GSH_Damp100.txt', 'Ala_fixed.txt', ...
    'Tau.txt', 'Ser.txt', 'Cys.txt'
};

%% Metabolic Composition
% Initialize default levels
% default_values = [... 
%     90,  30, ... % NAA+NAAG
%     40,  20, ... % Cr+PCr
%     18,   7, ... % GPC+PCh
%     25,  25, ... % Glu, Gln
%      5, 100, ... % Gly, Ins
%     30,  30, ... % GSH, Ala
%     0, 0, 0]; % Tau, Ser, Cys

% Metabolite.coefficients = [... 
%     0,  0, ...  % NAA+NAAG
%     10,  10, ...  % Cr+PCr
%     3,  4, ...  % GPC+PCh
%     10,  10, ...  % Glu, Gln
%     10, 10, ...  % Gly, Ins
%     10,  10, ...  % GSH, Ala
%     10,  10, 10 ]; % Tau, Ser, Cys

Metabolite.coefficients = [... 
    10, 2, ... % NAA+NAAG
    5,  1, ... % Cr+PCr
    2,  1, ... % GPC+PCh
    5,  2, ... % Glu, Gln
    5,  5, ... % Gly, Ins
    0,  0, ... % GSH, Ala
    0, 0, 0]; % Tau, Ser, Cys

% Initialize Metabolite.coefficient_structure using 'Metabolite.names' and 'Metabolite.coefficients'
Metabolite.coefficient_structure = cell2struct(num2cell(Metabolite.coefficients), Metabolite.names, 2);

% Variable Levels
Metabolite.variable_mets = {'Glu', 'Gln'};  % Metabolites to vary
Metabolite.variable_met_levels = {
    [5, 6]; 
    [3, 4];  
};

% % Variable Levels
% Metabolite.variable_mets = {'Ala', 'GSH'};  % Metabolites to vary
% Metabolite.variable_met_levels = {
%     [0, 10];  % Levels for 'Ala'
%     [0, 10];  % Levels for 'GSH'
% };

%% Tissue Types and Parameters
% noise levels. 
% 50: low noise, 55: okay, relatively low noise, 58: still okay
% 59: SNR 14, 60: complete mess
% filter levels
% 1: does almost nothing, % 2: looks okay
% 2.5: wrong fits with low CRLBs, but spectrum looks ok
% 3: completely broken
MetaInfo.noise_dbs = 0; % 59:0.25:60; %4e3:4e3:12.2e4;
MetaInfo.filter_levels = 1; %1.5:0.5:2.5;

%% CPU Cores
MetaInfo.CPU_cores = 4;

