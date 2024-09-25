%% FitCheck
% A simulation model for MRS data that aims to evaluate the fitting 
% accuracy of LCModel at different quality levels (SNR, FWHM, CRLB)

%% 1. Initialize
% Load default parameters
if exist('custom_parameters.m', 'file')
    disp('Loading custom parameters...');
    run('custom_parameters.m');
else
    disp('Loading default parameters...');
    run('default_parameters.m');
end

% Initialized dispstat
dispstat('','init'); 

reload_data = 0;

%% Add MATLAB functions to path (ExplonentialFilter, ...)
if ~contains(path, '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_MUSICAL_Streamlined_Git/MatlabFunctions'); ... 
        addpath(genpath('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_MUSICAL_Streamlined_Git/MatlabFunctions')); ...
end

%% 3. Load Data from Files
if reload_data == 1
    num_metabolites = length(met_files);
    metabolites = cell(1, num_metabolites);
    metabolite_signals_time = cell(1, num_metabolites);
    metabolite_signals_freq = cell(1, num_metabolites);

    dispstat(sprintf(''),'keepthis');
    dispstat(sprintf('Import basis spectra...'),'keepthis','timestamp');
    % Import data and process each metabolite
    for i = 1:num_metabolites
        dispstat(sprintf('%s', mets{i}));

        % Construct the file path
        filepath = fullfile(basis_sets_dir, met_files{i});

        % Import the data from the file
        data = importdata(filepath);

        % Store the metabolite name
        metabolites{i} = mets{i};

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


%% 4. Set Tissue Scaling Factors
% Metabolite coefficients are set in *_parameters.m

% Generate all combinations of variable levels using ndgrid
[level_grids{1:numel(mets_variable)}] = ndgrid(met_lvls_variable{:});

% Flatten the grids into vectors
for i = 1:numel(level_grids)
    level_grids{i} = level_grids{i}(:);
end

% Combine levels into combinations
lvls_combination = [level_grids{:}];
num_combinations = size(lvls_combination, 1);

% Preallocate tissue_scaling_factors with met_lvls_default fields
tissue_scaling_factors = repmat(met_lvls_default, num_combinations, 1);
tissue_comps = cell(num_combinations, 1);

% Assign variable levels to each combination and generate tissue_comps
for i = 1:num_combinations
    % Assign variable levels directly to the preallocated structure
    for m = 1:numel(mets_variable)
        metabolite = mets_variable{m};
        level_value = lvls_combination(i, m);
        tissue_scaling_factors(i).(metabolite) = level_value;
    end

    % Generate tissue type name by inspecting tissue_scaling_factors
    % We can create the name based on all metabolites that differ from default
    name_parts = {};
    for m = 1:numel(mets)
        metabolite = mets{m};
        default_value = met_lvls_default.(metabolite);
        current_value = tissue_scaling_factors(i).(metabolite);
        if ismember(metabolite,mets_variable) 
            name_parts{end+1} = sprintf('%s%g', metabolite, current_value);
        end
    end

    if isempty(name_parts) 
        tissue_comps{i} = 'Default';
    else
        tissue_comps{i} = strjoin(name_parts, '_');
    end
end

%% 5. Calculate parameters and preallocate InArray.csi
% Initialize data structure for storing results
num_noise_levels = length(noise_dbs);  
num_filter_widths = length(filter_levels);  

% Determine original, truncated and zerofilled signal lengths
siglen_orig = length(metabolite_signals_time{1});
siglen_trun = round(durmeas / dursim * siglen_orig);
siglen_new = 960;
siglen_zf = 2 * siglen_new;  

% InArray: Pre-allocate the field "csi" and add "DimNames"
InArray.csi = zeros(num_noise_levels, num_filter_widths, length(tissue_comps), siglen_zf);
InArray.DimNames = {'Noise', 'Linewidth', 'Composition', 'Spectrum'};

%% 6. Model signals
dispstat(sprintf('Creating CSI array.'),'keepthis','timestamp');

% Pick a tissue composition
for tissue_index = 1:length(tissue_comps)
    tissue_comp = tissue_comps{tissue_index};

    % Pick a noise levels
    noise_index = 0;
    for noise_db = noise_dbs
        noise_index = noise_index + 1;
       
        % Pick a filter widths
        filter_index = 0;        
        for filter_width_exp = filter_levels            
            filter_width = exp(filter_width_exp);
            filter_index = filter_index + 1;
            
            % Build up the signal as a linear combination of basis signals
            scalings = tissue_scaling_factors(tissue_index);
            signal = zeros(siglen_orig,1);
            for m = 1:num_metabolites
                signal = signal + metabolite_signals_time{m} * scalings.(metabolites{m});
            end
                       
            % Truncate and undersample signal
            signal = signal(1:siglen_trun);
            signal_us = interp1(1:siglen_trun, signal, linspace(1, siglen_trun, siglen_new));
            signal = signal_us;

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
%             fprintf('Noise dB: % 2i. Filter width: % 3i. Tissue: %s.\n', noise_db, filter_width, tissue_comp);

            dispstat(sprintf('Voxel: [%2i %2i %2i], Noise db: % 2.3e, Filter width: % 3.3e, Tissue: %s', ...
            noise_index, filter_index, tissue_index, noise_db, filter_width, tissue_comp), 'keepthis');
            
        end % linewidth loop
    end % noise loop
end % tissue loop

dispstat(sprintf('InArray.csi has been created!'), 'keepprev', 'timestamp');
fprintf('Dimensions: %i x %i x %i x %i\n', size(InArray.csi))

%% Debug

spectrum = real(fft(squeeze(InArray.csi(1,1,1,:))));
plot(spectrum);

% Define colors (3 colors with 3 shades each)
colors = {
    [1, 0.6, 0.6], [1, 0.3, 0.3], [1, 0, 0],    % Red shades
    [0.6, 1, 0.6], [0.3, 1, 0.3], [0, 1, 0],    % Green shades
    [0.6, 0.6, 1], [0.3, 0.3, 1], [0, 0, 1]     % Blue shades
};

figure; hold on;

% Loop over noise and filter levels
for i = 1 %:num_noise_levels
    for j = 1 %:num_filter_widths
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

