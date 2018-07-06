%% Clean up
clear all
close all
clc

%% Load data = [subject; run; trial; fix x; fix y; fix duration; stimulus; talker; static/dynamic; syllable (ba, ga, img)]
load('Clear_Syllables_Results.mat');
fix_report_corr_c = table2array(fix_report_corr);

%% Load ROI info
con_info = [1,1,1;2,1,2;3,1,3;4,2,1;5,2,2;6,2,3;7,3,1;8,3,2;9,3,3;10,4,1;11,4,2;12,4,3];
stim_info = {'MS1_F06_C1.mp4'; 'MS1_F06_C2.mp4'; 'MS1_F06_C3.mp4'; 'MS1_F08_C1.mp4'; 'MS1_F08_C2.mp4'; 'MS1_F08_C3.mp4'; 'MS1_M02_C1.mp4'; 'MS1_M02_C2.mp4'; 'MS1_M02_C3.mp4'; 'MS1_M03_C1.mp4'; 'MS1_M03_C2.mp4'; 'MS1_M03_C3.mp4';};
ROI_coords = load('ET_Dynamic_Static_McGurk_ROIs_07-24-2017_16-05-04.mat');
roi_list = {ROI_coords.ROI_coord.MS1_F06, ROI_coords.ROI_coord.MS1_F08,ROI_coords.ROI_coord.MS1_M02, ROI_coords.ROI_coord.MS1_M03};

%% Classify each fixation and measure distances to ROIs
check_rois = zeros(1,length(fix_report_corr_c));
dist_fix_roicenter = [];

display.dist = 90.0; % in cm
display.width = 100; % in cm
display.resolution = [ROI_coords.ROI_coord.image_monitor_details.monitor_res(1),ROI_coords.ROI_coord.image_monitor_details.monitor_res(2)]; % number of pixels of display in horizontal/vertical direction
pix_size = ROI_coords.ROI_coord.image_monitor_details.monitor_sizecm(1)/ROI_coords.ROI_coord.image_monitor_details.monitor_res(1); % cm/pix

for i = 1:length(fix_report_corr_c) % Run through all fixations
    
    % Select fix
    data = fix_report_corr_c(i,:); % info for 1 fixation
    talker = fix_report_corr_c(i,8);
    fix_x = data(4);
    fix_y = data(5);
        
    % Select ROI coords for talker
    roi = roi_list{talker}; % select roi values for this talker
    roi_vals = [roi.Eyes; roi.Mouth]; % get region structs for select roi as variables
    
    % Where is the fixation?
    for j = 1:numel(roi_vals) % 1 = Eyes; 2 = Mouth
        
        % Check if in ROI: 0 = not in roi; 1 = eyes; 2 = mouth
        if fix_x > roi_vals(j).pos_x(1) && fix_x < roi_vals(j).pos_x(2) && fix_y > roi_vals(j).pos_y(1) && fix_y < roi_vals(j).pos_y(2)
            check_rois(i) = j;
        end     
        
        % Check distance of mean fixation from ROI
        roi_center_x = roi_vals(j).center_x;
        roi_center_y = roi_vals(j).center_y;
        
        dist_roicenter_px = pdist([[fix_x fix_y];[roi_center_x roi_center_y]],'euclidean'); % Distance of mean fixation from ROI center
        dist_roicenter_cm = dist_roicenter_px*pix_size; % size in cm
        dist_roicenter_va = 2*180*atan(dist_roicenter_cm/(2*display(1).dist))/pi; % size in degree visual angle
        
        dist_fix_roicenter(i,j) = dist_roicenter_va;
        
    end
    
    % Check above/below line: 1 = upper; 0 = lower
    line_y = roi.Line.pos_y; 

    if fix_y <= line_y
        check_line(i) = 1; % upper
    else
        check_line(i) = 0; % lower
    end
        
end

roi_report_fix = [fix_report_corr_c check_rois' check_line' dist_fix_roicenter];

clearvars -except fix_report_corr roi_report_fix fix_report_corr_c trial_report_tot fix_report_tot

%% Summary per trial
sub = unique(roi_report_fix(:,1));

% CALCULATE: time in roi for each subject for each trial
roi_report_trial = [];
for i = 1:numel(sub)
    
    data_sub = roi_report_fix(roi_report_fix(:,1) == sub(i),:);
    
    runs = unique(data_sub(:,2));
    
    run_tot = [];
    for j = 1:numel(runs)
    
        data_run = data_sub(data_sub(:,2) == runs(j),:);     
        trials = unique(data_run(:,3));
        
        % Time in region during each trial
        trial_tot = [];
        for k = 1:numel(trials)
            data_trial = data_run(data_run(:,3) == trials(k),:);
            total_t = sum(data_trial(:,6));
            
            % Actual time in ms in each region
            data_eyes = sum(data_trial(data_trial(:,10) == 1,6));
            data_mouth = sum(data_trial(data_trial(:,10) == 2,6));
            data_out = sum(data_trial(data_trial(:,10) == 0,6));
            data_up = sum(data_trial(data_trial(:,11) == 1,6));
            data_low = sum(data_trial(data_trial(:,11) == 0,6));
            
            % Proportion of time in each region per trial
            prop_eyes = data_eyes/total_t;
            prop_mouth = data_mouth/total_t;
            prop_out = data_out/total_t;
            prop_low = data_low/total_t;
            prop_up = data_up/total_t;
            
            prop_tot = [prop_eyes prop_mouth prop_out prop_up prop_low];
            
            if size(data_trial,1) == 1
                trial_sum = data_trial(:,4:5);
                dist_sum = data_trial(:,12:13);
            else
                trial_sum = wmean(data_trial(:,4:5),data_trial(:,6));
                dist_sum = wmean(data_trial(:,12:13),data_trial(:,6));
            end

            num_fix = size(data_trial,1);
            trial_info = data_trial(1,1:3);
            con_info = data_trial(1,7:9);
            
            trial_tot(k,:) = [trial_info trial_sum total_t num_fix con_info prop_tot dist_sum];
        end
        
        run_tot = [run_tot;trial_tot];
    end

    roi_report_trial = [roi_report_trial;run_tot];
end

clearvars -except fix_report_corr roi_report_fix fix_report_corr_c trial_report_tot roi_report_trial fix_report_tot

%% Analysis by conditions
sub = unique(roi_report_fix(:,1));
data_sum = [];
for i = 1:numel(sub)
    
    data_s = roi_report_trial(roi_report_trial(:,1) == sub(i),:);
    data_s(:,13) = [];
    
    data_all = wmean(data_s(:,11:end),data_s(:,6));
    
    data_sum(i,:) = [sub(i) data_all];
    
end

stim = unique(roi_report_trial(:,8));
talk = unique(roi_report_trial(:,9));
syll = unique(roi_report_trial(:,10));
run = unique(roi_report_trial(:,2));

data_stim = [];
data_talk = [];
data_syll = [];
data_run = [];

for i = 1:numel(sub)
    
    data_s = roi_report_trial(roi_report_trial(:,1) == sub(i),:);
    
    for j = 1:numel(stim)
        stim_data = data_s(data_s(:,8)==j,:);
        data_stim(i,j) = nanmean(stim_data(:,15));        
    end
    
    for j = 1:numel(talk)
        talk_data = data_s(data_s(:,9)==j,:);
        data_talk(i,j) = nanmean(talk_data(:,15));        
    end
    
     for j = 1:numel(syll)
        syll_data = data_s(data_s(:,10)==j,:);
        data_syll(i,j) = nanmean(syll_data(:,15));        
     end

     for j = 1:numel(run)
        run_data = data_s(data_s(:,2)==j,:);
        data_run(i,j) = nanmean(run_data(:,15));        
     end
    
end

%% Save
head_line = {'Sub', 'Run', 'Trial', 'x', 'y', 'Fix_Duration', 'Stim', 'Talker', 'Syll', 'Region', 'Upper_lower', 'Dist_Eyes', 'Dist_Mouth'};
roi_report_fix = array2table(roi_report_fix, 'VariableNames', head_line);

head_line = {'Sub', 'Run', 'Trial', 'x', 'y', 'Fix_Duration', 'NumFix', 'Stim', 'Talker', 'Syll', 'PC_Eyes', 'PC_Mouth', 'PC_Off', 'PC_Up', 'PC_Low', 'Dist_Eyes', 'Dist_Mouth'};
roi_report_trial = array2table(roi_report_trial, 'VariableNames', head_line);

head_line = {'Sub' 'PC_Eyes', 'PC_Mouth', 'PC_Up', 'PC_Low', 'Dist_Eyes','Dist_Mouth'};
data_sum = array2table(data_sum, 'VariableNames', head_line);

save('ROI_Clear_Syllables_Results.mat','roi_report_fix','roi_report_trial','data_sum', 'data_stim', 'data_talk', 'data_syll','data_run');
