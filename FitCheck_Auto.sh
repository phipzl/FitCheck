#!/bin/bash

PARAMETERS_FILE="custom_parameters.m"

echo "Automatically generating custom parameters..."

# Check if custom_parameters.m already exists
if [ -f "$PARAMETERS_FILE" ] && [ -s "$PARAMETERS_FILE" ]; then
    echo "$PARAMETERS_FILE already exists and is not empty. Using existing file."
else
    echo "Creating $PARAMETERS_FILE..."
    cat <<EOL > $PARAMETERS_FILE
% custom_parameters.m

%% General Settings
reload_data = 1;
dispstat('', 'init');

%% Paths
% Note: These paths can be overridden by custom_config
lab_dir = '/your/lab/dir';
sim_dir = 'your/sim/dir';
process_results_dir = fullfile(lab_dir, 'your/process/results/dir');
output_dir = fullfile(lab_dir, 'your/output/dir');
batch_dir = fullfile(output_dir, 'BatchDir');
control_info_file = fullfile(sim_dir, 'your/control/info/file.m');
basis_sets_dir = fullfile(lab_dir, 'your/basis/sets/dir');

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
tissue_types = {'base'};
noise_dbs = 60;
filter_levels = 1.0;

%% Simulation Parameters
dursim = 1.092157;
durmeas = 0.3456;

%% MetaInfo Parameters
MetaInfo.DimNames = {'X'};
MetaInfo.Dimt1 = 4;
MetaInfo.pat_name = 'Meningioma_Simulation';
MetaInfo.LarmorFreq = 297.2232 * 1E6;
dwt = 3.2E-4;
MetaInfo.dwelltime = dwt * 1E9;

%% CPU Cores
CPU_cores = 1;
EOL
    echo "Custom parameters saved to $PARAMETERS_FILE."
fi

