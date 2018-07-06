clear all
close all
clc

%% Condition info
% NUMBER/SYLLABLE (1 = ba, 2 = da, 3 = ga, 4 = McG); TALKER
% con_info = [1 1; 2 2; 3 1; 4 1];
% stim_info = {'C1';'C2';'C3';'M1'};

% Prepare data matrices
center_report_tot = [];
center_report = [];
center_data_trial = [];

%% Check for previous results 
if exist('Clear_Syllables_Center_Fixes.mat','file') == 2 % -> already previous files
    
    check_data = 1;
    
    load('Clear_Syllables_Center_Fixes.mat');
    
    nb_sub_data = unique(table2array(center_report_tot(:,1))); % Number subjects already processed
    
    % Get all asc files
    result_files = dir('*Clear*.asc');
    result_files = {result_files(:).name}';
    
    % Get number of subjects not yet processed
    sub_files = zeros(length(result_files),1);
    
    for i = 1:length(result_files)
        
        sub_c = result_files{i};
        sub_files(i) = str2double(sub_c(1:2));
        
    end
    
    nb_sub_files = unique(sub_files);
    
    % Difference between data and files
    diff_files_data = setdiff(nb_sub_files,nb_sub_data);
    
    if isempty(diff_files_data) % in case all raw eye tracking data is already processed terminate script
        display('Nothing to do!');
        return
    end
    
    idx = find(sub_files >= min(diff_files_data));
    result_files = result_files(idx);
    
    %trial_report_tot_old = trial_report_tot;
    center_report_tot_old = table2array(center_report_tot);

else % not already previous files / first subject
    
    check_data = 0;
    result_files = dir('*Clear*.asc');
    result_files = {result_files(:).name}';
    
end

%% Read in & process data for each result file
for i = 1:numel(result_files) % For each result file that still needs to be processed (each sub & run is a separate file)
        
    % Make each line from file a string in it's own cell
    ET = strsplit(fileread(result_files{i}), '\n')'; 
    result_file_name = result_files{i};

    sub_num = str2double(result_file_name(1:2));
    idxr = strfind(result_file_name,'RUN');    
    run_num = str2double(result_file_name(idxr+3));

    fprintf('-- Subject %d, RUN %d --\n',sub_num,run_num)

    % Get sampling rate
    idx_sr = strfind(ET,'RECCFG');
    idx_sr = find(not(cellfun('isempty', idx_sr)));
    idx_sr = ET(idx_sr(1));
    idx_sr = cellfun(@(x) strsplit(x,' '),idx_sr(:),'uni',0);
    idx_sr = idx_sr{1,1};
    sampling_rate = str2double(cell2mat(idx_sr(4)));
    
    % Remove data from before SYNCTIME (expt start)
    idx = strfind(ET,'SYNCTIME');
    idx_sync = find(not(cellfun('isempty', idx)));
    ET = ET(idx_sync:end);
    
    % Get start and end for each center fix
    % Get indices for start & end of every movie
    idx = strfind(ET,'Center Fixation Start');
    idx_start = find(not(cellfun('isempty', idx)));
    idx = strfind(ET,'Center Fixation End');
    idx_end = find(not(cellfun('isempty', idx)));
    
    if sub_num == 1
        idx_start(1) = [];
        idx_end(1) = [];        
    end
    
    n_center = length(idx_start);
    
    for j = 1:n_center % run through data for each center cross
        
        fprintf('Center Fix Nr %d of %d\n',j,n_center)
        data_trial = ET(idx_start(j)+1:idx_end(j)-1); % Data for just this trial/stim
        
        % Clean up/remove first and last few samples
        data_trial(1:200) = []; % Remove first few samples since they are part of saccade prior to first fixation
        %data_trial(end-150:end) = []; % Remove last few samples since they are irrelevant

        % Get fixation start & end
        fix_start_idxv = strfind(data_trial,'SFIX'); % Starts of fixations
        fix_start_idx = find(not(cellfun('isempty', fix_start_idxv)));

        fix_end_idxv = strfind(data_trial,'EFIX'); % Ends of fixations
        fix_end_idx = find(not(cellfun('isempty', fix_end_idxv)));
        
        % CLEAN UP FIXATIONS
        % -> Outside movie
        if isempty(fix_start_idx) == 1 && numel(fix_end_idx) == 1
            fix_start_idx = 2;
        end
        if isempty(fix_end_idx) == 1 && numel(fix_start_idx) == 1
            fix_end_idx = length(data_trial)-1;
        end
        if isempty(fix_end_idx) == 1 && isempty(fix_start_idx) == 1
            fix_start_idx = 2;
            fix_end_idx = length(data_trial)-1;
        end
        
%         % -> If last fix is at very end of data
%         if fix_start_idx(end) > length(data_trial)-20
%             fix_start_idx(end) = [];
%         end
% 
%         % -> If first fix is at very start of data
%         if fix_end_idx(1) < 20
%            fix_start_idx(end) = [];
%         end
        
        if isempty(fix_start_idx) == 0 && isempty(fix_end_idx) == 0
            
        % -> Uneven number of starts & ends
            if numel(fix_start_idx) ~= numel(fix_end_idx)

        % ==> More starts than ends
                if numel(fix_start_idx) > numel(fix_end_idx)
            
                    if fix_start_idx(end) < round(length(data_trial)*0.95)
                        fix_end_idx = [fix_end_idx;length(data_trial)];
                    elseif fix_start_idx(end) >= round(length(data_trial)*0.95)
                        fix_start_idx(end) = [];
                    end
                    
        % ==> More ends than starts
                elseif numel(fix_start_idx) < numel(fix_end_idx)
                    
                    if fix_end_idx(1) > round(length(data_trial)*0.05)
                        fix_start_idx = [2; fix_start_idx];
                    elseif fix_end_idx(1) <= round(length(data_trial)*0.05)
                        fix_end_idx(1) = [];
                    end
                    
                end
        
        % -> Even number of starts & ends
             elseif numel(fix_start_idx) == numel(fix_end_idx)

        % ==> First end before first start
        % ---> Fixations at start AND end cut off
                if fix_end_idx(1) < fix_start_idx(1)
                        comp = [fix_start_idx fix_end_idx];

        % ----> Fixations at both ends long enough to be used 
                    if fix_start_idx(1) > round(length(data_trial)*0.05) && fix_end_idx(end) < round(length(data_trial)*0.95)
                        fix_start_idx = [2; fix_start_idx];
                        fix_end_idx = [fix_end_idx; length(data_trial)-1];
            
        % ----> Fixations at both ends too short = remove both
                    elseif fix_start_idx(1) <= round(length(data_trial)*0.05) && fix_end_idx(end) >= round(length(data_trial)*0.95)
                        fix_start_idx(end) = [];
                        fix_end_idx(1) = [];
                 
        % ----> Start ok; end too short
                    elseif fix_start_idx(1) > round(length(data_trial)*0.05) && fix_end_idx(end) >= round(length(data_trial)*0.95)
                        fix_start_idx(end) = [];
                        fix_start_idx = [2; fix_start_idx];
        
        % ----> Start too short; end ok
                    elseif fix_start_idx(1) <= round(length(data_trial)*0.05) && fix_end_idx(end) < round(length(data_trial)*0.95)
                        fix_end_idx(1) = [];
                        fix_end_idx = [fix_end_idx; length(data_trial)-1];
                    end
                end
            end
        
       else
           fix_start_idx = 2;
           fix_end_idx = length(data_trial) - 2;
        end
        
        center_data_trial = [];
        center_data_trial_mean = [];
        for p = 1:numel(fix_start_idx)
            
            fix_string = data_trial(fix_start_idx(p)+1:fix_end_idx(p)-1); % data for fix as string
            cells_fix = cellfun(@(x) strsplit(x,' '),fix_string(:),'uni',0); % split data into cells (by spaces)
            
            % Zeros matrix for each line x 2 (= number of values we want from each line)
            fix_pos = zeros(numel(cells_fix),2);
            
            for q = 1:numel(cells_fix) % loop through each line in fixation
                data_cur = cells_fix{q};
                pos_x = str2double(cell2mat(data_cur(2)));  % x position
                pos_y = str2double(cell2mat(data_cur(3)));  % y position
                %fix_dur(q,:) = str2double(cell2mat(data_cur(4)));
                fix_pos(q,:) = [pos_x pos_y];
            end
            
            % Take mean of lines to get fix pos (unless only 1)
            if size(fix_pos,1) == 1
                fix_mean = fix_pos;
            else
                fix_mean = nanmean(fix_pos);
                %dur_mean = nanmean(fix_dur);
            end
            
            fix_duration = size(fix_pos,1)/sampling_rate*1000;
            
            % Fixation means and durations in list
            center_data_trial(p,:) = [sub_num j fix_mean fix_duration];
        end
        
        if size(center_data_trial,1) > 1
            center_data_trial_mean = wmean(center_data_trial(:,3:4),center_data_trial(:,5));
        else
            center_data_trial_mean = center_data_trial(:,3:4);
        end
        
        center_data_trial_mean = [sub_num run_num j center_data_trial_mean];
        center_report = [center_report; center_data_trial_mean];
     
    end % End for this stimulus
    
end % for this result file

%% Save processed data
head_line = {'Sub' 'Run' 'Center_Fix' 'x' 'y'};

display('Save...');

if check_data == 1 % in case of add-on analysis
    center_report_tot = [center_report_tot_old; center_report];
    center_report_tot = array2table(center_report_tot,'VariableNames',head_line);
    save('Clear_Syllables_Center_Fixes.mat','center_report_tot');
else
    center_report_tot = array2table(center_report,'VariableNames',head_line);
    save('Clear_Syllables_Center_Fixes.mat','center_report_tot');
end

display('Analysis done!');
