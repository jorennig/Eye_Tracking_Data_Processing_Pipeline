clear all
close all
clc

load('Clear_Syllables_Results.mat'); % Fixation report 
fix_report = table2array(fix_report_tot);
load('Clear_Syllables_Center_Fixes.mat'); % Fixation report (during center crosses)
center_report = table2array(center_report_tot);

monitor_res = [1920,1080];
center_pos = monitor_res/2;

clearvars -except center_report fix_report fix_report_tot monitor_res center_pos monitor_res

corr_center = 0;

%% Calculate difference b/w measured central fix vs. "correct" central fix
sub = unique(fix_report(:,1));

if corr_center == 1

    run_diff = [];
    sub_diff = [];
    fix_report_corr = [];
    excl_sub = [];

    for i = 1:numel(sub) % Loop through each sub

        data_sub = fix_report(fix_report(:,1) == sub(i),:); % Get data for this sub (currently uncorrected)
        center_sub = center_report(center_report(:,1) == sub(i),:); % Get data for this sub
        run = unique(data_sub(:,2));

        for j = 1:numel(run)

            data_run = data_sub(data_sub(:,2) == run(j),:);
            center_run = center_sub(center_sub(:,2) == j,:);
            trials = unique(data_run(:,3));
            drift_idx = [15:15:numel(trials)];

            for k = 1:size(center_run,1) % run through each center cross

                % Get difference between subject's gaze pos during center fix &
                % actual position of center fix
                fix_c = center_run(center_run(:,3)==k,:);
                diff_cross = [center_pos - fix_c(:,4:5)]; % MONITOR center - SUBJECT gaze center

                % Use difference to correct trials before this fixation (and
                % after for last fixation
                if k == 1
                    data_corr = data_run(data_run(:,3)<=drift_idx(k),:);
                elseif k == size(center_run,1) && size(center_run,1) < numel(drift_idx)
                    data_corr1 = data_run(logical(double(data_run(:,3)>drift_idx(k-1)).*double(data_run(:,3)<=drift_idx(k))),:);
                    data_corr2 = data_run(data_run(:,3)>drift_idx(k),:);
                    data_corr = [data_corr1;data_corr2];
                else
                    data_corr = data_run(logical(double(data_run(:,3)>drift_idx(k-1)).*double(data_run(:,3)<=drift_idx(k))),:);
                end

                data_corr(:,4) = data_corr(:,4) + diff_cross(1);
                data_corr(:,5) = data_corr(:,5) + diff_cross(2);

                fix_report_corr = [fix_report_corr; data_corr];

            end
        end
    end

else
    
    fix_report_corr = fix_report;
    
end

%% Exclude fixations outside stimulus area
% Stimulus window parameters
win_x = [492 1428]; % x-vals (left,right)
win_y = [189 891]; % y-vals (low,high)

for i = 1:length(fix_report_corr) % loop through each subject
    
    data_fix = fix_report_corr(i,:);
    fix_x = data_fix(4);
    fix_y = data_fix(5);

    if fix_x > (win_x(2)) || fix_x < (win_x(1)) % Exclude if larger than high x val or smaller than low x val
        fix_report_corr(i,:) = NaN;
    end

    if fix_y > (win_y(2)) || fix_y < (win_y(1)) % outside of stimulus area
        fix_report_corr(i,:) = NaN;
    end

end

fix_report_corr(~any(~isnan(fix_report_corr), 2),:) = [];

%% Exclude fixations in off-center fixation areas
% Remove trials near off center fixation cross
fix_lh = struct('pos_x', [702 782],'pos_y', [299 379]); % left high
fix_ll = struct('pos_x', [702 782],'pos_y', [701 781]); % left low
fix_rh = struct('pos_x', [1113 1183],'pos_y', [299 379]); % right high
fix_rl = struct('pos_x', [1113 1183],'pos_y', [701 781]); % right low

fix_coords = [fix_lh; fix_ll; fix_rh; fix_rl];

for i = 1:numel(sub)
    
    data_sub = fix_report_corr(fix_report_corr(:,1) == sub(i),:);
    runs = unique(data_sub(:,2));
    
    for j = 1:numel(runs)
    
        data_run = data_sub(data_sub(:,2) == runs(j),:);     
        trials = unique(data_run(:,3));
        
        % Time in region during each trial
        for k = 1:numel(trials)
            data_trial = data_run(data_run(:,3) == trials(k),:);
            
            if size(data_trial,1) > 1
                
                fix_first = data_trial(1,:);
                fix_x = fix_first(4);
                fix_y = fix_first(5);
                
                for f = 1:size(fix_coords,1) % loop through each off center fix
                    
                    % Exclude any fixation on initial fixation crosses
                    if fix_x > fix_coords(f).pos_x(1) && fix_x < fix_coords(f).pos_x(2) && fix_y > fix_coords(f).pos_y(1) && fix_y < fix_coords(f).pos_y(2)
                        pos =  find(ismember(fix_report_corr,fix_first,'rows'));
                        fix_report_corr(pos,:) = NaN;
                    end
                    
                    % Exclude any fixation that are very short at the
                    % beginning of a trial
                    if fix_first(6) <= 100
                        pos =  find(ismember(fix_report_corr,fix_first,'rows'));
                        fix_report_corr(pos,:) = NaN;
                    end
                    
                end
                
            end
        end

    end 
    
end

fix_report_corr(~any(~isnan(fix_report_corr), 2),:)=[];

%% Save
head_line = {'Sub' 'Run' 'Trial' 'x' 'y' 'Fix_Duration' 'Stim' 'Talker' 'Syll'};
% Stim: stimulus numbers
% Talker: id talker (#1-5), 5 different talkers
% Syll: 1 = ba, 2 = ga, 3 = mcg/da

fix_report_corr = array2table(fix_report_corr,'VariableNames',head_line);

save('Clear_Syllables_Results.mat','fix_report_tot','fix_report_corr');
