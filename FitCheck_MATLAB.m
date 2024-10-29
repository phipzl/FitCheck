%% FitCheck
% A simulation model for MRS data that aims to evaluate the fitting 
% accuracy of LCModel at different quality levels (SNR, FWHM, CRLB)

%% 1. Initialize
% Initialized dispstat
dispstat('','init'); 

% Load default parameters
if exist('custom_parameters.m', 'file')
    dispstat(sprintf('Loading custom parameters...'),'timestamp');
    run('custom_parameters.m');
    dispstat(sprintf('Custom parameters loaded.'),'keepthis','timestamp');
else
    dispstat(sprintf('Loading default parameters...'),'timestamp');
    run('default_parameters.m');
    dispstat(sprintf('Default parameters loaded.'),'keepthis','timestamp');
end

% Add MATLAB functions to path (ExplonentialFilter, ...)
addpath(genpath('.'))
% if ~contains(path, '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_MUSICAL_Streamlined_Git/MatlabFunctions'); ... 
%         addpath(genpath('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Sourcecode/MRSI_Processing_ReleaseVersions/Part1_Reco_LCModel_MUSICAL_Streamlined_Git/MatlabFunctions')); ...
% end

%% 3. Load Data from Files
if MetaInfo.reload_data == 1 || ~isfield(Metabolite, 'import') || ~isfield(Metabolite.import, 'data')
    Metabolite.number = length(Metabolite.files);
    % Metabolites = cell(1, Metabolite.number);
    Metabolite.signals_time = cell(1, Metabolite.number);
    Metabolite.signals_freq = cell(1, Metabolite.number);

    % dispstat(sprintf(''),'keepthis');
    dispstat(sprintf('Import basis spectra...'),'timestamp');
    % Import data and process each Metabolite
    for i = 1:Metabolite.number
        dispstat(sprintf('%s', Metabolite.names{i}));

        Paths.filepath = fullfile(Paths.basis_sets_dir, Metabolite.files{i});
        Metabolite.import = importdata(Paths.filepath);

        Metabolite.signal_time = Metabolite.import.data(:,1) + 1j.*Metabolite.import.data(:,2);
        Metabolite.signals_time{i} = Metabolite.signal_time;

        % Perform FFT to convert to frequency domain
        Metabolite.signal_feq = fftshift(fft(Metabolite.signal_time));
        Metabolite.signals_freq{i} = Metabolite.signal_feq;
    end
    dispstat(sprintf('Basis spectra imported.'),'keepthis','timestamp');
else
    dispstat(sprintf('Basis spectra were not reloaded.'),'keepthis','timestamp');
end

% Optionally, convert cell arrays to matrices if all signals have the same length
% This is useful for easier manipulation and analysis
% time_signal_matrix = cell2mat(Metabolite.signals_time.');
% freq_signal_matrix = cell2mat(Metabolite.signals_freq.');


%% 4. Set Tissue Scaling Factors
% Simulation.metabolite coefficients are set in *_parameters.m
dispstat(sprintf('Setting tissue scaling factors...'),'timestamp');
clear Simulation % avoids mismatch issues when changing parameters for a second run

% If Metabolite.variable_mets is undefined, set the first Metabolite as variable and just use the default value
if ~exist('Metabolite.variable_mets', 'var')
    Metabolite.variable_mets = Metabolite.names(1);
    Metabolite.variable_met_levels = {Metabolite.coefficient_structure.(Metabolite.variable_mets{1})};
    dispstat(sprintf('Variable metabolites not defined. Using %s as variable.', Metabolite.variable_mets{1}),'keepthis','timestamp');
end

% Generate all combinations of variable levels using ndgrid
[Simulation.level_grids{1:numel(Metabolite.variable_mets)}] = ndgrid(Metabolite.variable_met_levels{:});

% Flatten the grids into vectors
for i = 1:numel(Simulation.level_grids)
    Simulation.level_grids{i} = Simulation.level_grids{i}(:);
end

% Combine levels into combinations
Simulation.lvls_combination = [Simulation.level_grids{:}];
Simulation.num_combinations = size(Simulation.lvls_combination, 1);

% Preallocate Simulation.tissue_scaling_factors with Metabolite.coefficient_structure fields
Simulation.tissue_scaling_factors = repmat(Metabolite.coefficient_structure, Simulation.num_combinations, 1);
Simulation.tissue_comps = cell(Simulation.num_combinations, 1);

% Assign variable levels to each combination and generate Simulation.tissue_comps
for i = 1:Simulation.num_combinations
    % Assign variable levels directly to the preallocated structure
    for m = 1:numel(Metabolite.variable_mets)
        Simulation.metabolite = Metabolite.variable_mets{m};
        Simulation.level_value = Simulation.lvls_combination(i, m);
        Simulation.tissue_scaling_factors(i).(Simulation.metabolite) = Simulation.level_value;
    end

    % Generate tissue type name by inspecting Simulation.tissue_scaling_factors
    % We can create the name based on all Simulation.metabolites that differ from default
    Simulation.name_parts = {};
    for m = 1:numel(Metabolite.names)
        Simulation.metabolite = Metabolite.names{m};
        Simulation.default_value = Metabolite.coefficient_structure.(Simulation.metabolite);
        Simulation.current_value = Simulation.tissue_scaling_factors(i).(Simulation.metabolite);
        if ismember(Simulation.metabolite,Metabolite.variable_mets) 
            Simulation.name_parts{end+1} = sprintf('%s%g', Simulation.metabolite, Simulation.current_value);
        end
    end

    if isempty(Simulation.name_parts) 
        Simulation.tissue_comps{i} = 'Default';
    else
        Simulation.tissue_comps{i} = strjoin(Simulation.name_parts, '_');
    end
end

dispstat(sprintf('Tissue scaling factors set.'),'keepthis','timestamp');

%% 5. Calculate parameters and preallocate InArray.csi
dispstat(sprintf('Creating CSI array.'),'keepthis','timestamp');

% Initialize data structure for storing results
Simulation.num_noise_levels = length(MetaInfo.noise_dbs);  
Simulation.num_filter_widths = length(MetaInfo.filter_levels);  

% Determine original, truncated and zerofilled signal lengths
Signal.OrigLength = length(Metabolite.signals_time{1});
Signal.TruncatedLength = round(MetaInfo.MeasDuration / MetaInfo.SimDuration * Signal.OrigLength);
Signal.NewLength = 1020;  % should be set by user
Signal.ZFLength = 2 *  Signal.NewLength;  

% InArray: Pre-allocate the field "csi" and add "DimNames"
InArray.csi = zeros(Simulation.num_noise_levels, Simulation.num_filter_widths, length(Simulation.tissue_comps), Signal.ZFLength);
InArray.DimNames = {'Noise', 'Linewidth', 'Composition', 'Spectrum'};

% create info.txt in Paths.out_dir
if exist(fullfile(Paths.out_dir, 'output.txt'), 'file')
    delete(fullfile(Paths.out_dir, 'output.txt'));
end
fileid = fopen(fullfile(Paths.out_dir, 'output.txt'), 'w');

% Write all parameters to the file
fprintf(fileid, 'Parameters:\n');
fprintf(fileid, 'OrigLength: %i\n', Signal.OrigLength);
fprintf(fileid, 'TruncatedLength: %i\n', Signal.TruncatedLength);
fprintf(fileid, 'NewLength: %i\n', Signal.NewLength);
fprintf(fileid, 'ZFLength: %i\n', Signal.ZFLength);
fprintf(fileid, 'Noise levels: %s\n', num2str(MetaInfo.noise_dbs));
fprintf(fileid, 'Filter levels: %s\n', num2str(MetaInfo.filter_levels));
fprintf(fileid, 'Tissue compositions: %s\n', strjoin(Simulation.tissue_comps, ', '));
fprintf(fileid, 'Metabolites: %s\n', strjoin(Metabolite.names, ', '));
% Write the Metabolite.coefficient_structure to file
fprintf(fileid, '--------------------------------\n');
fprintf(fileid, 'Metabolite base coefficients:\n');
for i = 1:numel(Metabolite.names)
    fprintf(fileid, '%s: %g\n', Metabolite.names{i}, Metabolite.coefficient_structure.(Metabolite.names{i}));
end
fprintf(fileid, '--------------------------------\n');
% Now for variable metabolites
    fprintf(fileid, 'Variable metabolites:\n');
for i = 1:numel(Metabolite.variable_mets)
    fprintf(fileid, '%s ', Metabolite.variable_mets{i});
    fprintf(fileid, '%g ', Metabolite.variable_met_levels{i});
    fprintf(fileid, '\n');
end
fprintf(fileid, '--------------------------------\n');
fprintf(fileid, '\n');

%% 6. Model signals

% Pick a tissue composition
for i_tissue = 1:length(Simulation.tissue_comps)
    Simulation.tissue_comp = Simulation.tissue_comps{i_tissue};

    % Pick a noise levels
    i_noise = 0;
    for noise_db = MetaInfo.noise_dbs
        i_noise = i_noise + 1;
       
        % Pick a filter widths
        filter_index = 0;        
        for filter_width_exp = MetaInfo.filter_levels            
            filter_width = exp(filter_width_exp);
            filter_index = filter_index + 1;
            
            % Build up the signal as a linear combination of basis signals
            scalings = Simulation.tissue_scaling_factors(i_tissue);
            signal = zeros(Signal.OrigLength,1);
            for m = 1:Metabolite.number
                signal = signal + Metabolite.signals_time{m} * scalings.(Metabolite.names{m});
            end
                       
            % Truncate and undersample signal 
            signal = signal(1:Signal.TruncatedLength);
            signal_us = interp1(1:Signal.TruncatedLength, signal, linspace(1, Signal.TruncatedLength,  Signal.NewLength));
            % Update:
            signal = signal_us;

            % Filter
            [signal, filter] = ExponentialFilter(signal, 360003, filter_width, 2);
            
            % Add White Gaussian Noise 
            %noise = sqrt(10^noise_db)/sqrt(2) * randn(1, 960) + 1i * sqrt(10^noise_db)/sqrt(2) * randn(1, 960); % because wgn is undefined (noise = wgn(1, 960, noise_db);)
            noise = wgn(1,Signal.NewLength,noise_db);
            %hold on; plot(noise)
            signal = signal + noise;

            % Zero-fill the signal
            signal_zf = zeros(1, Signal.ZFLength);
            signal_zf(1: Signal.NewLength) = signal(1: Signal.NewLength);

            % Flip spectrum by complex conjugating the time-domain signal
            signal_zf_conj = conj(signal_zf);

            % Store the result in the InArray structure
            InArray.csi(i_noise, filter_index, i_tissue, :) = signal_zf_conj;  
            
            % Output information  
%             fprintf('Voxel coordinates: % 3i % 3i % 3i. ', noise_index, filter_index, i_tissue);          
%             fprintf('Noise dB: % 2i. Filter width: % 3i. Tissue: %s.\n', noise_db, filter_width, Simulation.tissue_comp);

            dispstat(sprintf('Voxel: [%2i %2i %2i], Noise db: % 2.3e, Filter width: % 3.3e, Tissue: %s', ...
            i_noise, filter_index, i_tissue, noise_db, filter_width, Simulation.tissue_comp), 'keepthis');

            % Write the same info to a file
            fprintf(fileid, 'Voxel: [%2i %2i %2i], Noise db: % 2.3e, Filter width: % 3.3e, Tissue: %s\n', ...
            i_noise, filter_index, i_tissue, noise_db, filter_width, Simulation.tissue_comp);
            
        end % linewidth loop
    end % noise loop
end % tissue loop

% Close the file
fclose(fileid);

% Create mask for Write_LCM_files.
InArray.mask = ones(size(InArray.csi,1), size(InArray.csi,2), size(InArray.csi,3));   

dispstat(sprintf('InArray.csi has been created!'), 'keepprev', 'timestamp');
fprintf('Dimensions: %i x %i x %i x %i\n', size(InArray.csi))

%% Debug
% % Define colors for 9 spectra using the 'lines' colormap
% % colors = lines(size(InArray.csi, 3));  % Returns a 9x3 matrix of RGB values
% colors = [ 0 0 1 ; 0.3 0.7 1 ; 1 0 0; 1 0.7 0.3 ]; % Note: transposed!
% 
% ppm_scale = (compute_chemshift_vector(MetaInfo.LarmorFreq,... 
%     MetaInfo.dwelltime_seconds, size(InArray.csi, 4)));
% 
% % Loop over noise and filter levels
% for i = 1 %:Simulation.num_noise_levels
%     for j = 1 %:Simulation.num_filter_widths
%         figure; hold on;
%         for k = 1:size(InArray.csi, 3)
%             % Calculate the FFT and plot the magnitude of the FFT result
%             spectrum = abs(fftshift(fft(squeeze(InArray.csi(i,j,k,:)))));
%             plot(ppm_scale, spectrum, 'Color', colors(k,:), 'LineWidth', 1); 
%         end
% 
%         hold off;
%         set(gca, 'XDir', 'reverse');
%         xlabel('ppm'); ylabel('Magnitude');
%         title('Spectra of InArray.csi');
% 
%         % Use Simulation.tissue_comps directly as legend entries
%         legend(Simulation.tissue_comps, 'Interpreter', 'none');
%     end
% end

%% 7. Write LCM Files
close all

fprintf('##\n'); 
fprintf('Dwell time in us: \t% 10.2f \t(should be around 300 us)    \n', MetaInfo.dwelltime_seconds*1E6);
% fprintf('Bandwidth in Hz:  \t% 10.2f \t(should be around 3000 Hz)   \n', bandwidth);
% fprintf('Bandwidth in PPM: \t% 10.2f \t(should be around 10 ppm)    \n', 1E6 * bandwidth/MetaInfo.LarmorFreq);
fprintf('Orig length is: \t% 10i \t    \n', Signal.OrigLength);
fprintf('Truncated to: \t\t% 10i \t    \n', Signal.TruncatedLength);
fprintf('Downsampled to: \t% 10i \t    \n',  Signal.NewLength);
fprintf('Zero-filled to: \t% 10i \t    \n', Signal.ZFLength);
fprintf('##\nInArray.csi dims: % 2i % 2i % 2i % 2i \t   \n##\n', size(InArray.csi));




delete(fullfile(Paths.batchdir, '*'));
dispstat(sprintf('Writing LCModel files...'), 'keepprev', 'timestamp');
Write_LCM_files(InArray,Paths,MetaInfo,Paths.ControlInfo,InArray.mask,MetaInfo.CPU_cores);
dispstat(sprintf('Done.'), 'keepprev', 'timestamp');
fprintf('\n##################################################\n\n');

