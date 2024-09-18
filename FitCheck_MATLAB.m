%% A simulation model for MRS data that aims to evaluate the fitting accuracy of LCModel at different quality levels (SNR, FWHM, CRLB)


%% 1. Initialize
% Load default parameters
if exist('custom_parameters.m', 'file')
    disp('Loading custom parameters...');
    run('custom_parameters.m');
else
    disp('Loading default parameters...');
    run('default_parameters.m');
end


% Define paths
reload_data = 1;
lab_dir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab';
sim_dir = 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/Ala_GSH_Simulations'; 
process_results_dir = fullfile(lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/Results_v2/');
output_dir = fullfile(lab_dir, 'Process_Results/Tumor_Patients/Meningioma_Paper_2024/FittingSimulationResults');
batch_dir = fullfile(output_dir, 'BatchDir');
lcmodel_program_path = '/usr/local/lcmodel/bin/lcmodel';
%control_info_file = fullfile(sim_dir, 'ControlFiles/LCModel_Control_GH_2023_Meningeoma_new_v2.m');
control_info_file = fullfile(sim_dir, 'ControlFiles/LCModel_Control_GH_2020_Pat_sMM_ForRevision.m');
basis_sets_dir = fullfile(lab_dir, 'Basis_Sets/GH_FID_Basis_2019/jmrui');
basis_set =      fullfile(lab_dir, 'Basis_Sets/GH_FID_Basis_2019/Version_1_1_1/fid_0.000000ms.basis');
dispstat('','init'); % One time only initialization


%% Add MATLAB functions to path (ExplonentialFilter, ...)
if ~contains(path, '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_MUSICAL_Streamlined_Git/MatlabFunctions'); ... 
        addpath(genpath('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_MUSICAL_Streamlined_Git/MatlabFunctions')); ...
end

%% 2. Define metabolite files and names as well as tissue scaling factors
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

%% 3. Load Data from Files
if reload_data == 1
    num_metabolites = length(metabolite_files);
    metabolites = cell(1, num_metabolites);
    metabolite_signals_time = cell(1, num_metabolites);
    metabolite_signals_freq = cell(1, num_metabolites);

    dispstat(sprintf(''),'keepthis');
    dispstat(sprintf('Import basis spectra...'),'keepthis','timestamp');
    % Import data and process each metabolite
    for i = 1:num_metabolites
        dispstat(sprintf('%s', metabolite_names{i}));

        % Construct the file path
        filepath = fullfile(basis_sets_dir, metabolite_files{i});

        % Import the data from the file
        data = importdata(filepath);

        % Store the metabolite name
        metabolites{i} = metabolite_names{i};

        % Extract and store the time-domain signal
        signal_timedomain = data.data(:,1) + 1j.*data.data(:,2);
        metabolite_signals_time{i} = signal_timedomain;

        % Perform FFT to convert to frequency domain
        signal_freqdomain = fftshift(fft(signal_timedomain));
        metabolite_signals_freq{i} = signal_freqdomain;
    end
    dispstat(sprintf('Done!'),'keepthis','timestamp');
else
    fprintf('Data will not be loaded again!\n');
end

% Optionally, convert cell arrays to matrices if all signals have the same length
% This is useful for easier manipulation and analysis
% time_signal_matrix = cell2mat(metabolite_signals_time.');
% freq_signal_matrix = cell2mat(metabolite_signals_freq.');


%% 4. Set tissue scaling factors
% Tissue scaling factors with six cases (low/mid/high Ala and low/mid/high GSH)
tissue_scaling_factors = struct(...
    'base', struct('NAA', 1000, 'NAAG', 0, 'Cr', 0, 'PCr', 0, 'GPC', 0, 'PCh', 0, ...
                 'Glu', 0, 'Gln', 0, 'Gly', 0, 'Ins', 0, 'GSH', 0, 'Ala', 0, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'LL', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 10, 'GSH', 10, 'Ala', 10, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'LM', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 10, 'Ala', 30, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'LH', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 10, 'Ala', 50, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'ML', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 30, 'Ala', 10, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'MM', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 30, 'Ala', 30, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'MH', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 30, 'Ala', 50, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'HL', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 50, 'Ala', 10, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'HM', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 50, 'Ala', 30, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0), ...
    'HH', struct('NAA', 90, 'NAAG', 30, 'Cr', 40, 'PCr', 20, 'GPC', 18, 'PCh', 7, ...
                 'Glu', 25, 'Gln', 25, 'Gly', 5, 'Ins', 100, 'GSH', 50, 'Ala', 50, ...
                 'Tau', 0, 'Ser', 0, 'Cys', 0) ...
);



%% 5. Loop over different noise/linewidth parameters and different tissue types
% Define tissue types, noise and filter levels
tissue_types = {'base'}; %{'LL', 'LM', 'LH', 'ML', 'MM', 'MH', 'HL', 'HM', 'HH'};
noise_dbs = 60; % [ 10 20 30 ]; %4e3:4e3:12.2e4;
% noise levels
% 50: low noise
% 55: okay, relatively low noise
% 58: still okay
% 59: SNR 14
% 60: complete mess
filter_levels = 1.0; %[ 0.5 1.0 1.5 2.0] % [ 0.8 1.0 1.2 ]; %1.0:0.2:6; % 3 is too high, 1 does almost nothing. 
% filter levels
% 1: does almost nothing
% 2: looks okay
% 2.5: wrong fits with low CRLBs, but spectrum looks ok
% 3: completely broken


% Initialize data structure for storing results
num_noise_levels = length(noise_dbs);  
num_filter_widths = length(filter_levels);  

% Pre-calculate parameters for signal processing
dursim = 1.092157;     % simulation duration
durmeas = 0.3456;      % measurement duration

siglen_orig = length(metabolite_signals_time{1});
siglen_trun = round(durmeas / dursim * siglen_orig);
siglen_new = 960;
siglen_zf = 2 * siglen_new;  % Zero-filled length

% Pre-allocate the InArray.csi for efficiency
InArray.csi = zeros(num_noise_levels, num_filter_widths, length(tissue_types), siglen_zf);
%% 6. Model signals
dispstat(sprintf('Creating CSI array.'),'keepthis','timestamp');

% Loop over tissue types
for tissue_index = 1:length(tissue_types)
    tissue_type = tissue_types{tissue_index};
    % Loop over noise levels
    noise_index = 0;
    for noise_db = noise_dbs
        noise_index = noise_index + 1;
        
        
        % Loop over filter widths
        filter_index = 0;        
        for filter_width_in = filter_levels            
            filter_width = exp(filter_width_in);
            filter_index = filter_index + 1;
            
            % Create summed signals using scaling factors
            scalings = tissue_scaling_factors.(tissue_type);
            signal = zeros(size(metabolite_signals_time{1}));  % Initialize with the correct size            
            for m = 1:num_metabolites
                signal = signal + metabolite_signals_time{m} * scalings.(metabolites{m});
            end
                       
            % Truncate and undersample signal
            signal = signal(1:siglen_trun);
            signal_uns = interp1(1:siglen_trun, signal, linspace(1, siglen_trun, siglen_new));
            signal = signal_uns;

            % Filter
            [signal, filter] = ExponentialFilter(signal, 360003, filter_width, 2);
            
            % Add White Gaussian Noise 
            %noise = sqrt(10^noise_db)/sqrt(2) * randn(1, 960) + 1i * sqrt(10^noise_db)/sqrt(2) * randn(1, 960); % because wgn is undefined (noise = wgn(1, 960, noise_db);)
            noise =wgn(1,960,noise_db);
            %hold on; plot(noise)
            signal = signal + noise;

            % Zero-fill the signal
            signal_zf = zeros(1, siglen_zf);
            signal_zf(1:siglen_new) = signal(1:siglen_new);

            % Flip spectrum by complex conjugating the time-domain signal
            signal_zf_conj = conj(signal_zf);

            % Store the result in the InArray structure
            InArray.csi(noise_index, filter_index, tissue_index, :) = signal_zf_conj;  
            
            % Output information  
%             fprintf('Voxel coordinates: % 3i % 3i % 3i. ', noise_index, filter_index, tissue_index);          
%             fprintf('Noise dB: % 2i. Filter width: % 3i. Tissue: %s.\n', noise_db, filter_width, tissue_type);

            dispstat(sprintf('Voxel: [%2i %2i %2i], Noise db: % 2.3e, Filter width: % 3.3e, Tissue: %s', ...
            noise_index, filter_index, tissue_index, noise_db, filter_width, tissue_type), 'keepthis');
            
        end % linewidth loop
    end % noise loop
end % tissue loop

dispstat(sprintf('InArray.csi has been created!'), 'keepprev', 'timestamp');
fprintf('Dimensions: %i x %i x %i x %i\n', size(InArray.csi))

%% Debug

% Define colors (3 colors with 3 shades each)
colors = {
    [1, 0.6, 0.6], [1, 0.3, 0.3], [1, 0, 0],    % Red shades
    [0.6, 1, 0.6], [0.3, 1, 0.3], [0, 1, 0],    % Green shades
    [0.6, 0.6, 1], [0.3, 0.3, 1], [0, 0, 1]     % Blue shades
};

figure; hold on;

% Loop over noise and filter levels
for i = 1:num_noise_levels
    for j = 1:num_filter_widths
        % Compute the index for color selection
        color_index = (i-1)*length(filter_levels) + j;

        % Calculate the FFT and plot the magnitude of the FFT result
        spectrum = abs(fft(squeeze(InArray.csi(i,j,1,:))));
        plot(spectrum, 'Color', colors{color_index});
    end
end

hold off;
xlabel('Frequency Index'); ylabel('Magnitude');
legend('Red Shade 1', 'Red Shade 2', 'Red Shade 3', ...
       'Green Shade 1', 'Green Shade 2', 'Green Shade 3', ...
       'Blue Shade 1', 'Blue Shade 2', 'Blue Shade 3');
title('Spectra of InArray.csi Signals with Different Noise and Filter Levels');


%% 7. Write LCM Files
close all
MetaInfo.DimNames       = {'X'} % {'N','L','T'};    % noise, linewidth, tissue 
MetaInfo.Dimt1          = 4;                
MetaInfo.pat_name       = 'Meningioma_Simulation';

MetaInfo.LarmorFreq     = 297.2232 * 1E6;       % [LarmorFreq] = 1 Hz
% dwt                     = 1/2778;
dwt                     = 3.2E-4;
% dwt                     = 1/(297  * 10.5);      % Simulation: 100 ppm = 29722 Hz. -> dwt = 1/29722 Hz. 
                                                % Factor 10.8 from where? Basis set?                      

MetaInfo.dwelltime      = dwt * 1E9;            % [MetaInfo.dwelltime]: 1 ns for lcmodel
bandwidth               = 1/dwt;                % [dwt] = 1 s => [bandwidth] = 1 Hz

fprintf('##\n'); 
fprintf('Dwell time in us: \t% 10.2f \t(should be around 300 us)    \n', dwt*1E6);
fprintf('Bandwidth in Hz:  \t% 10.2f \t(should be around 3000 Hz)   \n', bandwidth);
fprintf('Bandwidth in PPM: \t% 10.2f \t(should be around 10 ppm)    \n', 1E6 * bandwidth/MetaInfo.LarmorFreq);
fprintf('Orig length is: \t% 10i \t    \n', siglen_orig);
fprintf('Truncated to: \t\t% 10i \t    \n', siglen_trun);
fprintf('Downsampled to: \t% 10i \t    \n', siglen_new);
fprintf('Zero-filled to: \t% 10i \t    \n', siglen_zf);
fprintf('##\nInArray.csi dims: % 2i % 2i % 2i % 2i \t   \n##\n', size(InArray.csi));


Paths.out_dir           = output_dir;
Paths.basis_file        = {basis_set};
Paths.LCM_ProgramPath   = lcmodel_program_path; 
Paths.batchdir          = batch_dir;

ControlInfo             = control_info_file;
mask                    = ones(size(InArray.csi,1), size(InArray.csi,2), size(InArray.csi,3));   
CPU_cores               = 1;

% Write_LCM_files(InArray,Paths,MetaInfo,ControlInfo,mask,CPU_cores)
% InArray:              Complex array with time domain data, e.g. [64x64x1x2048] CSI Matrix
% Paths:                Struct with Path info:
%                       -   Paths.out_dir 
%                       -   Paths.basis_file 
%                       -   Paths.LCM_ProgramPath 
%                       -   Paths.batchdir: Folder to which the bash scripts to start LCModel processing should be written
% MetaInfo:             Struct with metainfo about InArray: 
%                       -   MetaInfo.DimNames:       The dimension names of InArray. E.g. {'x','y','z'}. LCModel is instructed to create the output files as 
%                       "pat_name_DimNames1a_DimNames2b_DimNames3c" where a,b,c are numbers (e.g. the voxel [32 32 1]).
%                       -   MetaInfo.Dimt1:          Tells the function in which dimension of InArray the FID's are saved, i.e. which dimension is the
%                                                     t1-dimension / vecSize-dimension
%                       -   MetaInfo.pat_name:       The patient name. Used for the output file names
%                       -   MetaInfo.LarmorFreq:     The Larmor frequency which should be read out of the DICOM or raw-data header 
%                       -   MetaInfo.dwelltime:      The dwelltime, i.e. the duration between two consecutive time points.
% ControlInfo:          Struct or a string referring .m file with control parameters for LCModel fitting. This file / struct is executed at the end of the Control-Parameter settings,
%                       and can thus overwrite ALL control parameters. Be careful!
%                       You can also pass over a path with control parameters and additional parameters. The path is passed over by field .Path. The other parameters normal, e.g. .Others1 = ...
%                       In this case, specify which one should have priority by providing
%                       ControlInfo.Priority = 'Fields' or [...] = 'File'.
% mask:                 Only write files for voxels with 1 in mask; If all should be processed: mask = ones(...) (same size as InArray but no vecSize dim)
% CPU_cores:            Determines to how many batch scripts the individual voxels are split, so that those batch scripts can be run in parallel.

delete(fullfile(batch_dir, '*'));


dispstat(sprintf('Writing LCModel files...'), 'keepprev', 'timestamp');
Write_LCM_files(InArray,Paths,MetaInfo,ControlInfo,mask,CPU_cores);
dispstat(sprintf('Done.'), 'keepprev', 'timestamp');
fprintf('\n##################################################\n\n');

