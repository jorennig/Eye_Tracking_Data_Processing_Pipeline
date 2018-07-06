clear all
close all
clc

con_trials = [1,2,3];
mcg_trials = 4;

path_c = pwd;
path_plot = path_c;

load('/Volumes/data/BCM/Eye_Tracking/ROIs_Stimuli/ET_Noisy_ROIs_09-19-2017_12-21-07.mat');
monitor_res = ROI_coord.image_monitor_details.monitor_res;
center_position = monitor_res/2;
img_size_org = ROI_coord.image_monitor_details.image_size;

sf = ROI_coord.image_monitor_details.size_factor;
img_height = monitor_res(2)*sf;
sfac = img_height/img_size_org(2);
img_width = img_size_org(1)*sfac;

stim_size = [img_width img_height];

% % Check if there are .jpg files
% if ~isempty(dir('*.jpg'))
% 
%     check_data = 1;
%     
%     load('/Volumes/data/BCM/Eye_Tracking/ET_Data_Dynamic_Static/Dynamic_Static_Results.mat');
%     nb_sub_data = unique(fix_report_tot(:,1)); % Number subjects already processed
% 
%     % Get all asc files
%     data_files = dir('*Heatmap*.jpg');
%     for i=1:length(data_files)    
%         result_files{i} = data_files(i).name;       
%     end
%     clear data_files
% 
%     % Get number of subjects not yet processed
%     sub_files = zeros(length(result_files),1);
% 
%     for i = 1:length(result_files)
% 
%         sub_c = result_files{i};
%         idx = strfind(sub_c,'.');
%         sub_files(i) = str2double(sub_c(idx-2:idx-1));
% 
%     end
% 
%     % Difference between data and files
%     diff_files_data = setdiff(nb_sub_data,sub_files);
% 
%     if isempty(diff_files_data) % in case all raw eye tracking data is already processed terminate script
%         display('Nothing to do!');
%         return
%     end
%     
%     load('/Volumes/data/BCM/Eye_Tracking/ET_Data_Dynamic_Static/Dynamic_Static_Results.mat');
%     data_ds = table2array(fix_report_corr(ismember(fix_report_corr(:,1),diff_files_data),:));
%     load('/Volumes/data/BCM/Eye_Tracking/ET_Data_Dynamic_McGurk/McGurk_Dynamic_Results.mat');
%     data_dm = table2array(fix_report_corr(ismember(fix_report_corr(:,1),diff_files_data),:));
% 
% else
load('/Volumes/data/BCM/Eye_Tracking/ET_Data_Noisy/Noisy_Corr_Results.mat');
data = table2array(fix_report_corr);
% end

%% Prepare data
sub = unique(data(:,1));
data_re.data = data;

%% Create heat map
con = fieldnames(data_re);
title_stim = {'Noisy sentences'};

% Plot variables
lw = 1.1; % Line Width
bw = 0.8; % Barwidth
fsb = 16; % Font Size big
fsbb = 18; % Font Size small
path = path_plot; % Current path
rez = 300; % resolution (dpi) of final graphic
format = '-djpeg'; % image format

area = 20;

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))
    
    for j = 1:numel(con)
        
        % Heat map template
        heat_map = zeros(monitor_res(2),monitor_res(1));
        
        data_c = data_re.(con{j});        
        data_c = data_c(data_c(:,1)==sub(i),:);
        data_c = [round(data_c(:,4)) round(data_c(:,5)) data_c(:,6)];
        fix_overall = sum(data_c(:,3));
        
        data_c(data_c(:,1) > monitor_res(1),:) = NaN;
        data_c(data_c(:,2) > monitor_res(2),:) = NaN;
        data_c(any(isnan(data_c),2),:) = [];

        % Square area
        for k = 1:size(data_c,1)

            val = data_c(k,3);
            heat_map(data_c(k,2)-area:data_c(k,2)+area,data_c(k,1)-area:data_c(k,1)+area) = heat_map(data_c(k,2)-area:data_c(k,2)+area,data_c(k,1)-area:data_c(k,1)+area) + val;

        end
        
        % Calculate percent
        heat_map = heat_map/fix_overall*100;
        
        % Smooth heat map
        g_kernel = fspecial('gaussian',60,60);
        heat_map = filter2(g_kernel,heat_map);
        
        % Set threshold
        hm_values = sort(reshape(heat_map,[1 size(heat_map,1)*size(heat_map,2)]));
        hm_values(hm_values==0) = [];

        thr_pc = 90; % percent of lower values to be cut off
        vals = unique(hm_values);
        vals_sum = numel(vals);
        val_crit = round(vals_sum*(thr_pc/100));
        cut_off = vals(val_crit);

        % Apply threshold
        heat_map(heat_map < cut_off) = 0;

        heat_map = heat_map(monitor_res(2)/2-stim_size(2)/2:monitor_res(2)/2+stim_size(2)/2,monitor_res(1)/2-stim_size(1)/2:monitor_res(1)/2+stim_size(1)/2);
        %stim_size = ROI_coord.image_monitor_details.image_size;        
        %heat_map = heat_map(monitor_res(2)/2-stim_size(2)/2:monitor_res(2)/2+stim_size(2)/2-1,monitor_res(1)/2-stim_size(1)/2:monitor_res(1)/2+stim_size(1)/2-1);
        
%         % Plot    
%         cm_c = jet;
%         cm_c(1,:) = [1 1 1];
%         colormap(cm_c);
%         imagesc(heat_map);
% 
%         hold on
%         
%         img = rgb2gray(imread('MECOM1_001_pink_62dB.bmp'));
%         img = imresize(img,[img_height img_width]);
%         h = subimage(img);
%         set(h,'AlphaData',0.5);
%         %title(title_stim{j},'Fontsize',14);
% 
%         set(gca,'xtick',[]);
%         set(gca,'ytick',[]);
            
    end
    
%     % Save plot
%     set(gca,'LooseInset',get(gca,'TightInset'));
%     f = gcf; % f is the handle of the figure you want to export
%     if sub(i) < 10
%         name = ['Heatmap_Noisy_Sub0' num2str(sub(i))];
%     else
%         name = ['Heatmap_Noisy_Sub' num2str(sub(i))];
%     end
% 
%     set(gcf,'PaperPositionMode','auto')
%     print(f,fullfile(path,name),format,['-r',num2str(rez)])
% 
%     close all

end
