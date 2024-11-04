% default_parameters.m

%% General Settings
MetaInfo.reload_data = 0; % Set this to 1 to always reload the spectral data

%% Paths
Paths = struct();
Paths.lab_dir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab';
Paths.sim_dir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/FitCheck';
Paths.process_results_dir = fullfile(Paths.lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/Results_v2/');
Paths.out_dir = fullfile(Paths.lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/FittingSimulationResults');
Paths.batchdir = fullfile(Paths.out_dir, 'BatchDir');
Paths.LCM_ProgramPath = '/usr/local/lcmodel/bin/lcmodel';
%Paths.ControlInfo = fullfile(Paths.sim_dir, 'ControlFiles/LCModel_Control_GH_2020_Pat_sMM_ForRevision.m');
Paths.ControlInfo = fullfile(Paths.sim_dir, 'ControlFiles/LCModel_Control_PL_2024_ISMRM.m');
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
Metabolite.coefficients_WM = [... 
    25.1, 8.3, ... % NAA+NAAG
    0.0,  16.5, ... % Cr+PCr
    5.7,  0.0, ... % GPC+PCh
    14.2, 2.5, ... % Glu, Gln
    3.5,  19.9, ... % Gly, Ins
    3.3,  0, ... % GSH, Ala
    5.0, 3.6, 0]; % Tau, Ser, Cys

Metabolite.coefficients_GM = [... 
    29.3, 4.8, ... % NAA+NAAG
    0, 19.2, ... % Cr+PCr
    5.2, 0, ... % GPC+PCh
    27.7, 4.5, ... % Glu, Gln
    2.6, 19.7, ... % Gly, Ins
    1.7, 0, ... % GSH, Ala
    5.2, 2.8, 0]; % Tau, Ser, Cys


Metabolite.coefficients = Metabolite.coefficients_WM;

% Initialize Metabolite.coefficient_structure using 'Metabolite.names' and 'Metabolite.coefficients'
Metabolite.coefficient_structure = cell2struct(num2cell(Metabolite.coefficients), Metabolite.names, 2);

% Variable Levels
Metabolite.variable_mets = {'Gln','Gly'};  % Metabolites to vary
Metabolite.variable_met_levels = {
    [2.5 3.5 4.5];
    [3.5 4.5 5.5]
};

% % Variable Levels
% Metabolite.variable_mets = {'Ala', 'GSH'};  % Metabolites to vary
% Metabolite.variable_met_levels = {
%     [0, 10];  % Levels for 'Ala'
%     [0, 10];  % Levels for 'GSH'
% };

%% Parameters
MetaInfo.noise_dbs = 15:1:30; 
MetaInfo.filter_levels = 3.5:0.5:4.5; 



