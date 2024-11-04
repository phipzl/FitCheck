% function [MetData_relevant_amp MetData_relevant_sd] = Evaluate_Simulation_Results(out_dir, matrix_size)
% MetMaps is a 4D array (noise values x line width x tissue type x metabolite)

out_dir =['/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Process_Results/Tumor_Patients/Meningioma_Paper_2024/FittingSimulationResults'];
file_prefix = 'My_Simulation';
% matrix_size_override = [30 6];

%% Extracting matrix_size from file name: 
% List all .table files in the directory and sort them in descending order
table_files = dir(fullfile(out_dir, '*.table'));

[~, idx] = sort({table_files.name});  % Sort in ascending order
sorted_files = table_files(flip(idx));  % Flip to get descending order

% Extract numbers from the first (largest) filename
file_name = sorted_files(1).name;
tokens = regexp(file_name, 'X(\d+)_Y(\d+)_Z(\d+)', 'tokens');

if ~isempty(tokens)
    matrix_size = [str2double(tokens{1}{1}), str2double(tokens{1}{2}), str2double(tokens{1}{3})];
else
    error('Could not extract X, Y, Z values from the filename.');
end

if exist('matrix_size_override','var')
    matrix_size = matrix_size_override;
end

mask = ones(matrix_size);     


%% Create list of tables

cur_dir = cd(sprintf('%s',out_dir));                                    % save current directory and then go to directory out_dir_spectra
table_list = struct2cell(dir(sprintf('%s/*.table',out_dir)));           % dir(folder) creates a struct of all matching files; the structure fields are: name, date, bytes, isdir, datenum; create a cell out of struct
table_list = table_list(1,:);                                                   % only the first entry (the former .name, that is now the first row) is of interest.
cd(sprintf('%s',cur_dir));                                                      % go back to cur_dir directory

    

%% MetData_amp_title etc.
no_met = 0;
MetData_amp_title = {};
table_list_index = 0;
while(no_met == 0)          % if metabolites are found --> stop, otherwise open next table file

    table_list_index = table_list_index + 1;
    voxel_raw_in = sprintf('%s/%s', out_dir, table_list{table_list_index});     % pick an element of the list
    fid_vox = fopen(voxel_raw_in,'r');                                          % open this table-file
    
    % Read the number of metabolites
    tline = fgetl(fid_vox);
    while(~feof(fid_vox))         
      	if(strfind(tline, '$$CONC'))
            met_str = regexp(tline,' ','split');
            met_str(cellfun(@isempty,met_str)) = [];
            no_met = str2double(met_str{2}) - 1;            
%            if(WaterRefState == 3)	
%                no_met = str2double(met_str{2}) - 7;
%            else
%                no_met = str2double(met_str{2}) - 1;
%            end

            % read name of each metabolite line by line
            %%%%%%%%% fixed bug: if LCM gives output like "0.1+Scyllo" then the met_name could not be read out %%%%%%%%%%%%%%%
            fgetl(fid_vox);
            for n = 1:no_met
                tline = fgetl(fid_vox);

                met_line = regexp(tline,'\s+','split');                 % necessary because there can be a very small or large number in met_line{4} displayed as E-10;
                met_line(4)=strrep(met_line(4),'E-','yyx');             % in next step met_line(4) gets split by + or - deliminator but it should not be split when E- or E+;      
                met_line(4)=strrep(met_line(4),'E+','yyy');             % so E- and E+ gets replaced by xx and xy; after splitting, xx and xy get converted back to E- and E+ resp.

                if (regexp(met_line{4}, '+')>=0)                        % if met_line{4} contains a +, eg "0.328+Scyllo" (which should be met_line{4}=0.328 and met_line{5}=Scyllo) 
                    met_line(4:5) = regexp(met_line{4}, '+', 'split');  % split the string met_line{4} by + and write them into met_line(4) and (5)
                    met_line=[met_line, {''}];                          % add a sixth cell array element containing nothing, just to get the same cell array as in the case when the if doesn't match.
                elseif (regexp(met_line{4}, '-')>=0)                    % same if met_line{4} contains -; both (met_line{4} containts + and -) can never happen.
                    met_line(4:5) = regexp(met_line{4}, '-', 'split');  
                    met_line=[met_line, {''}];
                end

                met_line(4)=strrep(met_line(4),'yyx','E-');
                met_line(4)=strrep(met_line(4),'yyy','E+');

                MetData_amp_title = [MetData_amp_title; met_line{5}];
            end %for n      
        end % if strfind 
        
        
        tline = fgetl(fid_vox);
    end
    fclose(fid_vox);         % close the table file
    
end

if(no_met == 0)      % if no table with some metabolites was found, quit
    display([char(10) 'No table with any metabolite inside was found. The PC will explode now. Leave work and drink a tea.' char(10)])   %char(10) = new line
    quit force
end
            
MetData_sd_title = MetData_amp_title;
MetData_extra_title = [ {'FWHM'}; {'SNR'}; {'shift'}; {'0_pha'}; {'1_pha'}];
%% 5. READ METABOLITE DATA OF ALL TABLES

%preallocate memory for data
MetData_amp = NaN(size(mask,1),size(mask,2),size(mask,3),numel(MetData_amp_title));
MetData_sd = NaN(size(mask,1),size(mask,2),size(mask,3),numel(MetData_sd_title));
% MetData_amp_clipped = NaN(size(mask,1),size(mask,2),size(mask,3),numel(MetData_amp_title));
% MetData_sd_clipped = NaN(size(mask,1),size(mask,2),size(mask,3),numel(MetData_sd_title));
MetData_extra = NaN(size(mask,1),size(mask,2),size(mask,3),numel(MetData_extra_title));

fprintf('Read data...\n');

%read table files
voxel_no = 0;
for T=1:size(mask,3)            % tissue, 'z'
    for L=1:size(mask,2)        % line width, 'y'
        for N=1:size(mask,1)    % noise level, 'x'
            voxel_raw_in = sprintf('%s/%s_X%02d_Y%02d_Z%02d.table', out_dir, file_prefix, N, L, T);
            
            voxel_no = voxel_no+1;
            
            fid_vox = fopen(voxel_raw_in,'r');
            
            if (fid_vox == -1)
                
                MetData_amp(N,L,T,:) = 0;
                MetData_sd(N,L,T,:) = 0;				
                MetData_extra(N,L,T,:) = 0;
                continue;
                
            else
                
                tline = fgetl(fid_vox);
                while ischar(tline)
                    % find out number of metabolites
                    if(strfind(tline, '$$CONC'))
                        met_str = regexp(tline,' ','split');
                        met_str(cellfun(@isempty,met_str)) = [];
                        no_met = str2double(met_str{2}) - 1;
%                         if(WaterRefState == 3)
%                             no_met = str2double(met_str{2}) - 7;
%                         else
%                             no_met = str2double(met_str{2}) - 1;
%                         end

                        tline = fgetl(fid_vox);                    

                        % read data for each metabolite line by line
                        for n = 1:no_met
                            tline = fgetl(fid_vox);
                            met_line = regexp(tline,'\s+','split');
                            MetData_amp(N,L,T,n) = str2double(met_line{2});
                            MetData_sd(N,L,T,n)  = str2double(regexprep(met_line{3},'%',''));
                        end
                    end

                    % find out number of lines with extra info
                    if(strfind(tline, '$$MISC'))
                        misc_str = regexp(tline,' ','split');
                        misc_str(cellfun(@isempty,misc_str)) = [];
                        no_misc = str2num(misc_str{2}) - 1;

                        % read data for extra info line by line
                        tline = fgetl(fid_vox);
                        misc_line = regexp(tline,'\s+','split');

                        MetData_extra(N,L,T,1) = str2double(misc_line{4});    %FWHM
                        MetData_extra(N,L,T,2) = str2double(misc_line{8});    %SNR

                        tline = fgetl(fid_vox);
                        misc_line = regexp(tline,'=','split');
                        misc_line = regexp(misc_line{2},'ppm','split');

                        MetData_extra(N,L,T,3) = str2double(misc_line{1});    %shift

                        tline = fgetl(fid_vox);
                        misc_line = regexp(tline,':','split');
                        misc_line = regexp(misc_line{2},'\s+','split');

                        if isnan(str2double(misc_line{1}))
                            MetData_extra(N,L,T,4) = str2double(misc_line{2});    %zero-order phase
                            MetData_extra(N,L,T,5) = str2double(misc_line{4});    %first-order phase
                        else
                            MetData_extra(N,L,T,4) = str2double(misc_line{1});    %zero-order phase
                            MetData_extra(N,L,T,5) = str2double(misc_line{3});    %first-order phase
                        end
                    end

                    tline = fgetl(fid_vox);
                end

                fclose(fid_vox);
            end
        end
    end
end


%solves problem with NaN in cases where
% 1.no fitting possible
% 2.too large number !!!
% MetData_amp(isnan(MetData_amp))=0;
% MetData_sd(isnan(MetData_sd))=0;
 
% clearvars -except MetData_amp MetData_amp_title MetData_sd MetData_sd_title MetData_extra MetData_extra_title out_dir ...
% no_met mask mask_res CRLB_treshold_value LogMat_FWHM LogMat_SNR LogMat_combined MetData_amp_test vecSize print_individual_spectra_flag ...
% LarmorFreq dwelltime

%% Extract and organize data per tissue type and metabolite
% Initialize cell arrays for each metabolite
% MetData_norm = cell(4, 1); % Normalized data for each metabolite
% MetData_abs = cell(4, 1);  % Absolute data for each metabolite
% 
% for met = 1:4
%     MetData_norm{met} = squeeze(MetData_norm_amp(:, :, :, met));
%     MetData_abs{met} = squeeze(MetData_abs_amp(:, :, :, met));
% end

%% Normalize to reference
% What is the reference? -> noise_db of ~40, filter width of ~300
% This corresponds to: N = 4, L = 4.
Nref=1; Lref=1;
% Prepare a matrix of reference values for each tissue type and metabolite.
% Care: If there is only one tissue type, we would get a Nx1 matrix instead 
% of a 1x4 matrix, so we need to transpose in that case.
if size(MetData_amp,3) == 1 
    reference_matrix=squeeze(MetData_amp(Nref,Lref,:,:))';
else
    reference_matrix=squeeze(MetData_amp(Nref,Lref,:,:)); 
end

% Matrix: Tissue type (TT) x metabolite (Met)
% Tensor: SNR x LW x TT x Met 
reference_tensor = repmat(reference_matrix, [1 1 matrix_size(1) matrix_size(2)]);
reference_tensor = permute(reference_tensor, [3 4 1 2]);            

% Explanation: We use repmat because for the SNR and LW because these
% should use the same reference values! TT and Mets should have different 
% reference values!

% Normalize amplitudes
MetData_RelativeChange = (MetData_amp - reference_tensor) ./ reference_tensor;
% MetData_abs_amp = abs(MetData_amp - reference_tensor);


%% Define most interessting metabolites

target_metabolites = ["Glu", "Gln", "Gly", "Ins"];
[~, met_indices] = ismember(target_metabolites, string(MetData_amp_title));

if any(met_indices == 0)
    error('Some target metabolites were not found in MetData_amp_title.');
end



%% Extract SNR and line width (LW) values for each tissue type
LW = squeeze(MetData_extra(:, :, :, 1));
SNR = squeeze(MetData_extra(:, :, :, 2));

%% Check: Comparing Cases. Mets: Glu (r), Gln (g), Ins (b), Gly (c). 
% 
% noisel = 1;
% linewi = 1;
% 
% figure; hold on;
% plot(squeeze(MetData_amp(:, linewi, 1, 1)), 'r', 'DisplayName', MetData_title{1});
% plot(squeeze(MetData_amp(:, linewi, 1, 2)), 'g', 'DisplayName', MetData_title{2});
% plot(squeeze(MetData_amp(:, linewi, 1, 3)), 'b', 'DisplayName', MetData_title{3});
% plot(squeeze(MetData_amp(:, linewi, 1, 4)), 'c', 'DisplayName', MetData_title{4});
% 
% xlabel('Noise Level');
% ylabel('Amplitude');
% title(sprintf('Metabolite Fit Stability Against Noise (LW: %d)', linewi));
% 
% legend show; % Display the legend with metabolite titles and colors
% hold off;


%% Visualization
close all

% Select tissue type for plotting (e.g., z = 1)
tissue_idx = 1;

% Extract the corresponding SNR and LW values from MetData_extra
SNR_vals = MetData_extra(:,:,tissue_idx, 2);  % 3D tensor for SNR
LW_vals = MetData_extra(:,:,tissue_idx, 1);   % 3D tensor for LW (FWHM)

% Flatten SNR and LW values
SNR_flat = SNR_vals(:);
LW_flat = LW_vals(:);

% Calculate and display SNR and LW statistics (independent of metabolite)
SNR_stats = [min(SNR_flat), mean(SNR_flat), median(SNR_flat), max(SNR_flat), numel(unique(SNR_flat))];
LW_stats = [min(LW_flat), mean(LW_flat), median(LW_flat), max(LW_flat), numel(unique(LW_flat))];

fprintf('');
fprintf('SNR: Min = %2.3f, Mean = %2.3f, Median = %2.3f, Max = %2.3f, Unique values = %d\n', SNR_stats);
fprintf('LW:  Min = %2.3f, Mean = %2.3f, Median = %2.3f, Max = %2.3f, Unique values = %d\n', LW_stats);

% Loop through each metabolite for plotting
n=0;     
figure('Name','Metabolic Fit vs SNR and LW', ...
    'units','normalized','outerposition',[0 0 1 1])
for metabolite_idx = met_indices 
    n=n+1;
    
    % Extract the relevant 2D slice of MetData_amp for the selected tissue and metabolite
    data_slice = MetData_RelativeChange(:,:,tissue_idx,metabolite_idx);
    data_flat = data_slice(:);
    
    % Create a scatter plot with color-coded MetData_amp values
    subplot(2,2,n);
    scatter(LW_flat, SNR_flat, 60, data_flat, 'filled'); % 60 is the point size
    
    colorbar;
    caxis([min(data_flat), max(data_flat)]);
    colormap(jet);
    clim([-1 1]);  % Colorbar limits
    xlabel('Line Width (LW)');
    ylabel('Signal-to-Noise Ratio (SNR)');
    title(sprintf('Stability of Metabolic Fit for %s', MetData_amp_title{metabolite_idx}));
    grid on;
    set(gca, 'FontSize', 12);

end

set(gcf, 'Position', get(0, 'Screensize'));


%% Write Met Data to Files

% fprintf('Write to files...\n');
% 
% mkdir(out_dir, 'evaluation/');
% 
% for met = 1:size(MetData_amp,4)
%     % Orig Amplitudes (*_amp_map.raw)
%     path = sprintf('%s/evaluation/%s_amp_map.raw',out_dir,MetData_amp_title{met,1});
%     fid = fopen(path,'w+');
%     fwrite(fid,MetData_amp(:,:,:,met),'float');
%     fclose(fid);
% end

%% Raw To Minc (CURRENTLY BUGGED!)

% % Define the evaluation folder and retrieve list of .raw files
% evaluation_folder = fullfile(out_dir, 'evaluation');
% raw_files = dir(fullfile(evaluation_folder, '*.raw'));
% 
% % Extract dimensions dynamically from matrix_size
% dimensions_str = sprintf('%d ', matrix_size);
% dimensions_str = strtrim(dimensions_str);  % Remove any trailing space
% 
% % Loop through each .raw file and convert to .mnc
% for i = 1:length(raw_files)
%     raw_filename = fullfile(evaluation_folder, raw_files(i).name);
%     output_filename = fullfile(evaluation_folder, [raw_files(i).name(1:end-4), '.mnc']);
% 
%     % Build and execute the rawtominc command
%     command = sprintf('rawtominc -clobber -transverse -input "%s" "%s" %s', raw_filename, output_filename, dimensions_str);
%     system(command);
% end
% 
% fprintf('Done.\n');
% 
% 
% %end

