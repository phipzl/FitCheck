% default_parameters.m

%% General Settings
reload_data = 1;
dispstat('', 'init');

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

%% Metabolite Names and Files
metabolite_names = {
    'NAA', 'NAAG', 'Cr', 'PCr', 'GPC', 'PCh', ...
    'Glu', 'Gln', 'Gly', 'Ins', 'GSH', 'Ala', ...
    'Tau', 'Ser', 'Cys'
};

metabolite_files = {
    'NAA_damp55.txt', 'NAAG.txt', 'Cr_Damp190.txt', 'PCr_Damp190.txt', 'GPC.txt', 'PCh.txt', ...
    'glu.txt', 'Gln_Damp130.txt', 'glycine.txt', 'Ins.txt', 'GSH_Damp100.txt', 'Ala_fixed.txt', ...
    'Tau.txt', 'Ser.txt', 'Cys.txt'
};

%% Tissue Types and Parameters
tissue_types = {'base'};  % You can add more tissue types as needed
noise_dbs = 60;
filter_levels = 1.0;

%% Simulation Parameters
dursim = 1.092157;     % Simulation duration
durmeas = 0.3456;      % Measurement duration

%% MetaInfo Parameters
MetaInfo.DimNames = {'X'};
MetaInfo.Dimt1 = 4;
MetaInfo.pat_name = 'Meningioma_Simulation';
MetaInfo.LarmorFreq = 297.2232 * 1E6;  % [Hz]
dwt = 3.2E-4;                          % Dwell time in seconds
MetaInfo.dwelltime = dwt * 1E9;        % [ns] for LCModel

%% CPU Cores
CPU_cores = 1;

