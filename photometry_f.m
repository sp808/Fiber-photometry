%% Function declaration

% To-do list:
% - Remove plot displays mid-run (i.e., just save the plots to a file)
% - Work out a master Excel sheet or Matlab file that can store organized across runs

function [dFF, dFF_HC, prepost_dFF_discrete, summary_Mean, summary_Bouts, summary_prepost] = photometry_f(BLOCKPATH_input,BLOCKPATH_HC_input,onset_FF_input,filename_input,sheetname_input)

% Manually import session information

% Use this section if not running in batch function mode, comment out otherwise
close all; clear all; clc;

BLOCKPATH = 'vhpc_t_PFC_NAc_social_bar-180202-135210/195-181206-155238'
BLOCKPATH_HC = 'vhpc_t_PFC_NAc_social_bar-180202-135210/195-181206-154500'
onset_FF = 6
sheetname = '195_24 hour'
filename = 'Photometry #2'

%% Import data and define manual inputs, comment out if not running in function mode
% 
% BLOCKPATH = char(BLOCKPATH_input)
% BLOCKPATH_HC = char(BLOCKPATH_HC_input)
% onset_FF = str2double(onset_FF_input) %Manually input
% 
% filename = char(filename_input) % Photometry #2 - 193_21 day - 271_30 day - Event Logs // Photometry #1 - 269_24 hour - 577_30 day - Event Logs
% sheetname = char(sheetname_input)

%% Now read the specified data from our block into a Matlab structure.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'});

%% Import, normalize, trim, and synchronize CFC recording,

%create time vectors for each data stream by dividing by the sampling frequency
time_g22 = (1:length(data.streams.g2_2.data))/data.streams.g2_2.fs; 
time_uv22 = (1:length(data.streams.uv22.data))/data.streams.uv22.fs;

% Linear fit to isobestic + subtraction

bls = polyfit(data.streams.uv22.data(1:end), data.streams.g2_2.data(1:end), 1); %return the coefficients for a polynomial p(x) of degree n that is the best fit for the data in y

% so y = mx + b, where bls(1) = m, b = bls(2) (bls is just a vector that stores the polynomial coefficients in descending order
% so basically this is creating a vector of all 465 values fitted to this line? or something?
Y_fit_all = bls(1) .* data.streams.uv22.data + bls(2); 
Y_dF_all = data.streams.g2_2.data - Y_fit_all; %subtract element wise from the raw 465 signal

% dF/F

dFF = Y_dF_all ./ Y_fit_all;
line = zeros(size(Y_dF_all));

% Syncronize timestamp to FreezeFrame light pulse and trim

offset_FF = onset_FF + 300; %start time is manually input, stop time is 5 min later

% [min_value_onset, start_idx] = min(abs(time_g22 - onset_FF)); %subtract element-wise by the timestamped value and then find the index of the lowest value (closest to zero)
% [min_value_offset, stop_idx] = min(abs(time_g22 - offset_FF));

start_idx = ceil(onset_FF * data.streams.g2_2.fs);
stop_idx = ceil(offset_FF * data.streams.g2_2.fs);

dFF_trim = NaN(1, round(300 * data.streams.g2_2.fs));

k = 1;
for i = start_idx:stop_idx
    dFF_trim(k) = dFF(i);
    k = k + 1;
end

%% Import freezing times from Excel and determine freezing indices

read_Excel = xlsread(filename, sheetname,'H:H'); %Reads the relative time column from the Observer output Excel sheet
[m,n] = size(read_Excel);
num_Bouts = m / 2; %Assuming each start event precedes a stop event, the size of the this array / 2 is the number of bouts

start_times = NaN(num_Bouts, 1);
stop_times = NaN(num_Bouts, 1);
start_indices = NaN(num_Bouts, 1);
stop_indices = NaN(num_Bouts, 1);

f = 1;
ff = 1;

% first, split the start/stop time events into two parallel arrays (so an
% indice corresponds to pairs of start/stops in each array)
for i = 1:m
    
    ans = mod(i,2);
    
    if ans == 1
        
        start_times(f) = read_Excel(i,1);
        f = f + 1;
        
    elseif ans == 0
        
        stop_times(ff) = read_Excel(i,1);
        ff = ff + 1;
    end
end

% Add the absolute time of the CFC recording start to align to Synapse
start_times = start_times + onset_FF;
stop_times = stop_times + onset_FF;

% Next, convert times to index arrays in the original time vector
for i = 1:num_Bouts
    [m1,n1] = min(abs(time_g22 - start_times(i))); %subtract element-wise by the timestamped value and then find the index of the lowest value (closest to zero)
    [m2,n2] = min(abs(time_g22 - stop_times(i)));
    
    start_indices(i) = n1;
    stop_indices(i) = n2;
    
end

% Now, fill an actual array with photometry signal values only within
% freezing bouts

within_Freezing_dFF = NaN(size(dFF)); % photometry signal within freezing bouts
within_Freezing_dFF_discrete = NaN(num_Bouts, size(dFF,2)); % array to hold the photometry signal of each discrete bout in separate rows for subsequent analysis
within_Freezing_dFF_discrete_sorted = NaN(num_Bouts, size(dFF,2)); %array to hold bouts > 2 seconds (and other criteria)
outside_Freezing_dFF = NaN(size(dFF)); % photometry signal outside of bouts
Freezing_bouts = zeros(size(dFF)); % freezing bouts binary for computations
Freezing_visual = NaN(size(dFF)); %freezing bouts NaN vs. 0 for plotting
Freezing_visual2 = NaN(size(dFF)); %fill negative y values in plot

%fill dFF values within freezing bouts
for i = 1:num_Bouts
    for j = start_indices(i):stop_indices(i)
        within_Freezing_dFF(j) = dFF(j);
    end
end

within_Freezing_dFF((stop_idx+1):length(within_Freezing_dFF)) = []; %trim beginning and end of recording  
within_Freezing_dFF(1:(start_idx-1)) = []; 

%fill dFF values within freezing bouts, but with each bout being in a discrete row 
k1 = 1;
k2 = 1;

for i = 1:num_Bouts
    for j = start_indices(i):stop_indices(i)
        within_Freezing_dFF_discrete(i,k1) = dFF(j);
      
        if ((stop_indices(i) - start_indices(i)) / data.streams.g2_2.fs) >= 2
            
            within_Freezing_dFF_discrete_sorted(i,k2) = dFF(j);
            k2 = k2 + 1;    
        end
        k1 = k1 + 1;
    end
end

%create binary array for freezing behavior
for i = 1:num_Bouts
    for j = start_indices(i):stop_indices(i)
        Freezing_bouts(j) = 1;
    end
end

%visualize freezing bouts for plots
for i = 1:num_Bouts
    for j = start_indices(i):stop_indices(i)
        Freezing_visual(j) = max(dFF_trim); %can set this to zero, I am setting it to just above the peak dFF
        Freezing_visual2(j) = min(dFF_trim);
    end
end

Freezing_visual((stop_idx+1):length(Freezing_visual)) = []; %trim beginning and end of recording  
Freezing_visual(1:(start_idx-1)) = []; 
Freezing_visual2((stop_idx+1):length(Freezing_visual2)) = []; %trim beginning and end of recording  
Freezing_visual2(1:(start_idx-1)) = []; 

%fill dFF values for non-freezing bouts, doing this by filling in indices
%where Freezing_bouts = 0 (aka no freezing)

for i = start_idx:stop_idx
    if Freezing_bouts(i) == 0
        outside_Freezing_dFF(i) = dFF(i);
    end
end

outside_Freezing_dFF_untrimmed = outside_Freezing_dFF; % saving a copy of the untrimmed array since it gets used for analysis later
outside_Freezing_dFF((stop_idx+1):length(outside_Freezing_dFF)) = []; %trim beginning and end of recording  
outside_Freezing_dFF(1:(start_idx-1)) = [];  

%% Import HC recording, normalize, trim, and plot

data_HC = TDTbin2mat(BLOCKPATH_HC, 'TYPE', {'epocs', 'scalars', 'streams'});

time_g22_HC = (1:length(data_HC.streams.g2_2.data))/data_HC.streams.g2_2.fs; 
time_uv22_HC = (1:length(data_HC.streams.uv11.data))/data_HC.streams.uv22.fs;

% Linear fit to isobestic + subtraction (HC)

bls_HC = polyfit(data_HC.streams.uv22.data(1:end), data_HC.streams.g2_2.data(1:end), 1); %return the coefficients for a polynomial p(x) of degree n that is the best fit for the data in y

% so y = mx + b, where bls(1) = m, b = bls(2) (bls is just a vector that stores the polynomial coefficients in descending order
% so basically this is creating a vector of all 465 values fitted to this line? or something?
Y_fit_all_HC = bls_HC(1) .* data_HC.streams.uv22.data + bls_HC(2); 
Y_dF_all_HC = data_HC.streams.g2_2.data - Y_fit_all_HC; %subtract element wise from the raw 465 signal

% dF/F (HC)

dFF_HC = Y_dF_all_HC ./ Y_fit_all_HC;
line2 = zeros(size(Y_dF_all_HC));

% Trim middle 5 minutes of dF/F

dFF_trim_HC = NaN(size(dFF_HC));
[m,n] = (size(dFF_HC));
middle_idx = n./2;
start_idx_HC = middle_idx - 100000;
stop_idx_HC = middle_idx + 100000;

for i = start_idx_HC:stop_idx_HC
    dFF_trim_HC(i) = dFF_HC(i);
end

%% Extract and store across session CFC. HC parameters

summary_Mean(1) = nanmean(dFF_trim);
summary_Mean(2) = nanmin(dFF_trim);
summary_Mean(3) = nanmax(dFF_trim);
summary_Mean(4) = nanstd(dFF_trim);

summary_Mean(5) = nanmean(dFF_trim_HC);
summary_Mean(6) = nanmin(dFF_trim_HC);
summary_Mean(7) = nanmax(dFF_trim_HC);
summary_Mean(8) = nanstd(dFF_trim_HC);

summary_Mean(9) = nanmean(within_Freezing_dFF);
summary_Mean(10) = nanmin(within_Freezing_dFF);
summary_Mean(11) = nanmax(within_Freezing_dFF);
summary_Mean(12) = nanstd(within_Freezing_dFF);

summary_Mean(13) = nanmean(outside_Freezing_dFF);
summary_Mean(14) = nanmin(outside_Freezing_dFF);
summary_Mean(15) = nanmax(outside_Freezing_dFF);
summary_Mean(16) = nanstd(outside_Freezing_dFF);

%% Find peaks/transients

% tshold = (summary_Mean(3) - summary_Mean(2)) * 0.3; % min peak height based on Kim et al., 2019, 30% of the max - min height
% 
% % width threshold of 0.01 sec is arbitrary, consider changing in the future
% 
% [p] = findpeaks(dFF_trim, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01);  %find total number of peaks to make a new array
% peaks_dFF = NaN(4,length(p)); 
% %parse height, location, width, and prominence of peaks
% [peaks_dFF(1,:), peaks_dFF(2,:), peaks_dFF(3,:), peaks_dFF(4,:)] = findpeaks(dFF_trim, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01);
% 
% %do the same thing for HC 
% tshold = (summary_Mean(6) - summary_Mean(5)) * 0.3;
% 
% [p] = findpeaks(dFF_trim_HC, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01); 
% peaks_dFF_HC = NaN(4,length(p)); 
% [peaks_dFF_HC(1,:), peaks_dFF_HC(2,:), peaks_dFF_HC(3,:), peaks_dFF_HC(4,:)] = findpeaks(dFF_trim_HC, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01);
% 
% %do the same thing for within freezing
% tshold = (summary_Mean(8) - summary_Mean(9)) * 0.3;
% 
% [p] = findpeaks(within_Freezing_dFF, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01); 
% peaks_Frz = NaN(4,length(p)); 
% [peaks_Frz(1,:), peaks_Frz(2,:), peaks_Frz(3,:), peaks_Frz(4,:)] = findpeaks(within_Freezing_dFF, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01);
% 
% %do the same thing for mobile segment
% tshold = (summary_Mean(11) - summary_Mean(12)) * 0.3;
% 
% [p] = findpeaks(outside_Freezing_dFF, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01); 
% peaks_Mob = NaN(4,length(p)); 
% [peaks_Mob(1,:), peaks_Mob(2,:), peaks_Mob(3,:), peaks_Mob(4,:)] = findpeaks(outside_Freezing_dFF, data.streams.g2_2.fs, 'MinPeakHeight', tshold, 'MinPeakWidth', 0.01);
% 
% %calculate summary data for peaks, number of peaks, mean/max * height/width for CFC/HC/Frz/Mob (4 rows x 5 columns)
% summary_Peaks = NaN(1,20);
% summary_Peaks(1) = size(peaks_dFF,2);
% summary_Peaks(2) = nanmean(peaks_dFF(1,:));
% summary_Peaks(3) = nanmean(peaks_dFF(3,:));
% summary_Peaks(4) = max(peaks_dFF(1,:));
% summary_Peaks(5) = max(peaks_dFF(3,:));
% 
% summary_Peaks(6) = size(peaks_dFF_HC,2);
% summary_Peaks(7) = nanmean(peaks_dFF_HC(1,:));
% summary_Peaks(8) = nanmean(peaks_dFF_HC(3,:));
% summary_Peaks(9) = max(peaks_dFF_HC(1,:));
% summary_Peaks(10) = max(peaks_dFF_HC(3,:));
% 
% summary_Peaks(11) = size(peaks_Frz,2);
% summary_Peaks(12) = nanmean(peaks_Frz(1,:));
% summary_Peaks(13) = nanmean(peaks_Frz(3,:));
% summary_Peaks(14) = max(peaks_Frz(1,:));
% summary_Peaks(15) = max(peaks_Frz(3,:));
% 
% summary_Peaks(16) = size(peaks_Mob,2);
% summary_Peaks(17) = nanmean(peaks_Mob(1,:));
% summary_Peaks(18) = nanmean(peaks_Mob(3,:));
% summary_Peaks(19) = max(peaks_Mob(1,:));
% summary_Peaks(20) = max(peaks_Mob(3,:));

% generate transient landmarks for plots

% plot_peaks = NaN(1,size(dFF_trim,2));
% 
% for i = 1:size(peaks_dFF,2)
%     j1 = round(peaks_dFF(2,i) * data.streams.g2_2.fs);
%     j2 = round(peaks_dFF(2,i) * data.streams.g2_2.fs) + round(peaks_dFF(3,i) * data.streams.g2_2.fs);
%     
%         for k = j1:j2
%             plot_peaks(k) = 0; %dFF_trim (k);
%         end
% end

% for i = 1:size(peaks_dFF,2)
%     plot_peaks(1,round(peaks_dFF(2,i) * data.streams.g2_2.fs)) = 0;
% end

%% Plots

figure(1);

time_vector_300s = linspace(0,300,size(dFF_trim,2));
line3 = zeros(1, size(dFF_trim,2));

ax1 = subplot(4,2,1);
plot(time_g22, data.streams.g2_2.data(1,:),'Color', [0, 0.4470, 0.7410]);
hold on
plot(time_uv22, data.streams.uv22.data(1,:), 'k');
hold off
%axis tight;
%title({'GCaMP/UV raw signal (CFC)'}, 'FontSize',10);
ylabel('mV');
xlabel('Time (seconds)');
ax1.FontSize = 20;
%xlabel('Time (seconds)');

ax2 = subplot(4,2,3);
plot(time_g22, Y_dF_all(1,:),'k');
%title('Fitted GCaMP signal (CFC)', 'FontSize',10);
ylabel('dF');
%xlabel('Time (seconds)');

ax3 = subplot(4,2,5);
plot(time_g22, dFF(1,:),'k');
hold on
plot(time_g22, line,'k:');
hold off
%title('dF/F (CFC)', 'FontSize',10);
%ylim([-0.1 0.2])
ylabel('dF/F');
%xlabel('Time (seconds)');

ax4 = subplot(4,2,7);
plot(time_vector_300s, dFF_trim(1,:),'k');
hold on
plot(time_vector_300s, line3,'k:');
%plot(time_vector_300s, plot_peaks,'r', 'LineWidth', 2);
hold off
%title('dF/F, trimmed (CFC)', 'FontSize',10);
%ylim([-0.1 0.2])
ylabel('dF/F');
xlabel('Time (seconds)');

ax5 = subplot(4,2,2);
plot(time_g22_HC, data_HC.streams.g2_2.data(1,:), 'Color', [0, 0.4470, 0.7410]);
hold on
plot(time_uv22_HC, data_HC.streams.uv22.data(1,:),'k');
hold off
axis tight;
%title({'GCaMP/UV raw signal (HC)'}, 'FontSize',10);
ylabel('mV');
%xlabel('Time (seconds)');

ax6 = subplot(4,2,4);
plot(time_g22_HC, Y_dF_all_HC(1,:),'k');
%title('Fitted GCaMP signal (HC)', 'FontSize',10);
ylabel('dF');
%xlabel('Time (seconds)');

ax7 = subplot(4,2,6);
plot(time_g22_HC, dFF_HC(1,:),'k');
ylim([-0.1 0.2])
hold on 
plot(time_g22_HC, line2,'k:');
hold off
%title('dF/F (HC)', 'FontSize',10);
ylabel('dF/F');
%xlabel('Time (seconds)');

ax8 = subplot(4,2,8);
plot(time_g22_HC, dFF_trim_HC(1,:),'k');
hold on
plot(time_g22_HC, line2,'k:');
hold off
%title('dF/F, trimmed (HC)', 'FontSize',10);
ylim([-0.1 0.2])
ylabel('dF/F (%)');
xlabel('Time (seconds)');

figure(2);

Freezing_visual2 = Freezing_visual2 + 0.03;

axis1 = subplot(3,1,1);
plot(time_vector_300s, dFF_trim(1,:),'k');
hold on
plot(time_vector_300s, Freezing_visual2, 'LineWidth', 10, 'Color', [0, 0.4470, 0.7410]);
%area(time_vector_300s, Freezing_visual, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'ColorSpec', [0, 0.4470, 0.7410]);
%area(time_vector_300s, Freezing_visual2, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'ColorSpec', [0, 0.4470, 0.7410]);
plot(time_vector_300s, line3,'k:');
hold off
%title('dF/F, trimmed (CFC)', 'FontSize',10);
ylabel('dF/F (%)');
ylim([min(dFF_trim)-0.05 max(dFF_trim)+0.05])
%xlabel('Time (seconds)');

axis2 = subplot(3,1,2);
plot(time_vector_300s, within_Freezing_dFF(1,:),'k');
hold on
plot(time_vector_300s, line3,'k:');
hold off
%title('dF/F, freezing bouts', 'FontSize',10);
ylabel('dF/F (%)');
%xlabel('Time (seconds)');

axis3 = subplot(3,1,3);
plot(time_vector_300s, outside_Freezing_dFF(1,:),'k');
hold on
plot(time_vector_300s, line3,'k:');
hold off
%title('dF/F, non-freezing', 'FontSize',10);
ylabel('dF/F (%)');
xlabel('Time (seconds)');

%% Plots, per bout

% Bout-by-bout heatmap

max_boutLength = 0;

%find maximum bout length 
for i = 1:num_Bouts
    k = 1;
    for j = 1:size(dFF,2)
        if isnan(within_Freezing_dFF_discrete_sorted(i,j)) == 0
            k = k + 1;
        end
    end
    if max_boutLength < k
        max_boutLength = k;
    end
end

% align all sorted bouts to t = 0

within_Freezing_dFF_discrete_sorted_noNaN = NaN(num_Bouts, max_boutLength);

for i = 1:num_Bouts
    k = 1;
    for j = 1:size(dFF,2)
        if isnan(within_Freezing_dFF_discrete_sorted(i,j)) == 0
            within_Freezing_dFF_discrete_sorted_noNaN(i,k) = within_Freezing_dFF_discrete_sorted(i,j);            
            k = k + 1;
        end
    end
end

within_Freezing_dFF_discrete_sorted_noNaN(~any(~isnan(within_Freezing_dFF_discrete_sorted_noNaN), 2), :) = []; %remove rows containing only NaN (aka filtered bouts <2 sec)

% Create a column vector of size max_boutLength but only with values corresponding to every 5 seconds
x_axis_Labels = NaN(max_boutLength, 1);
y_axis_Labels = NaN(size(within_Freezing_dFF_discrete_sorted_noNaN,1),1);

k = 0;
    for i = 1:max_boutLength
        if round(i/data.streams.g2_2.fs) == k
        x_axis_Labels(i,1) = k;
        k = k + 5;
        end
    end

figure(3);

axes1 = subplot(3,1,1);

h = heatmap(within_Freezing_dFF_discrete_sorted_noNaN);
%h.Title = 'dF/F, individual bouts';
%h.XLabel = 'Time (seconds)';
h.YLabel = 'Bout';
h.Colormap = parula;
h.GridVisible = 'off';
h.ColorMethod = 'none';
h.MissingDataLabel = 'Non-freezing';
h.XDisplayLabels = x_axis_Labels;
h.YDisplayLabels = y_axis_Labels;

% Averaged bout trace with SEM shading

% create another matrix with the dFF of sorted bouts but only for the first 2 seconds
within_Freezing_dFF_discrete_sorted_2sec = NaN(size(within_Freezing_dFF_discrete_sorted_noNaN, 1), round(2 * data.streams.g2_2.fs));
time_vector_2sec = linspace(0, 2, round(2 * data.streams.g2_2.fs)); 

% fill array

for i = 1:size(within_Freezing_dFF_discrete_sorted_noNaN, 1)
    for j = 1:round(2 * data.streams.g2_2.fs)
        within_Freezing_dFF_discrete_sorted_2sec(i,j) = within_Freezing_dFF_discrete_sorted_noNaN(i,j);
    end
end

% heatmap of first two seconds

x_axis_Labels = NaN(size(time_vector_2sec,2), 1);
k = 0;

    for i = 1:size(time_vector_2sec,2)
        if round((time_vector_2sec(1,i)),1) == k
            x_axis_Labels(i,1) = k;
            k = k + 0.5;
        end
    end

axes2 = subplot(3,1,2);

h2 = heatmap(within_Freezing_dFF_discrete_sorted_2sec);
%h2.Title = 'dF/F over time, per bout (first 2 seconds)';
%h.XLabel = 'Time (seconds)';
h2.YLabel = 'Bout';
h2.Colormap = parula;
h2.GridVisible = 'off';
h2.ColorMethod = 'none';
h2.MissingDataLabel = 'Non-freezing';
h2.XDisplayLabels = x_axis_Labels;
h2.YDisplayLabels = y_axis_Labels;

line3 = zeros(size(within_Freezing_dFF_discrete_sorted_2sec, 2), 1);

axes3 = subplot(3,1,3);
stdshade(within_Freezing_dFF_discrete_sorted_2sec, 0.10, 'k', time_vector_2sec);
hold on 
plot(time_vector_2sec, line3, 'k:');
hold off
%title('Averaged dF/F across bouts, first 2 seconds', 'FontSize',10);
ylabel('dF/F');
xlabel('Time (seconds)');

%% Pre/post (2 seconds) bout analysis

% First task, of the > 2 second freezing bouts already isolated, find those that have at least 2 sec of mobility preceding it. 
% So we need to go back to within_Freezing_dFF which preserves the overall time structure, but we lose individual bout identity.
% So probably the approach is to do the original sorting method into discrete rows but with the added criterion of 2 sec of preceding mobility (which we have info on
% via within outside_Freezing_dFF

prepost_dFF = NaN(1, size(dFF,2));
prepost_dFF_discrete = NaN(num_Bouts, round(4 * data.streams.g2_2.fs));

for i = 1:num_Bouts
k3 = 1;
        if ((stop_indices(i) - start_indices(i)) / data.streams.g2_2.fs) >= 2 % first, check that the bout length is > 2 seconds, same as before
            
            % next, iterate through outside_Freezing_dFF with the start/stop indices offset 2 seconds backwards and check that none of these values are NaN (i.e., represent 
            % full bouts of mobility) and store this information in a boolean array. Then, do a final check of the boolean array itself and then fill in the dFF values
            % 2 seconds pre + post if everything checks out.
            
            mobility_Check = false(1, size(dFF,2)); % create/reset boolean array

            for j = (start_indices(i) - round(2 * data.streams.g2_2.fs)) : (start_indices(i))
                if outside_Freezing_dFF_untrimmed(1,j) == NaN   
                    mobility_Check(1,j) = 1;
                end
            end
            
            if any(mobility_Check) == 0
               for w = (start_indices(i) - round(2 * data.streams.g2_2.fs)) : (start_indices(i) + round(2 * data.streams.g2_2.fs))
                   prepost_dFF(1,w) = dFF(1,w);
                   prepost_dFF_discrete(i, k3) = dFF(1, w); %fill in discrete rows at the same time
                   k3 = k3 + 1;
               end
            end
        end
end

prepost_dFF((stop_idx+1):length(prepost_dFF)) = []; %trim beginning and end of recording  
prepost_dFF(1:(start_idx-1)) = [];  
prepost_dFF_discrete(any(isnan(prepost_dFF_discrete), 2), :) = []; %trim NaN-only rows

figure(4);

axe1 = subplot(3,1,1);
plot(time_vector_300s, prepost_dFF(1,:),'k');
hold on
plot(time_vector_300s, Freezing_visual, 'LineWidth', 10);
hold off
%title('pre/post dFF pairs', 'FontSize',10);
ylabel('dF/F (%)');
%xlabel('Time (seconds)');

% Pre/post heatmap of bout transition

time_vector_4sec = linspace(-2, 2, size(prepost_dFF_discrete,2)); 
x_axis_Labels = NaN(size(time_vector_4sec,2), 1);
y_axis_Labels = NaN(size(prepost_dFF_discrete,1),1);

k = -2;
    for i = 1:size(time_vector_4sec,2)
        if round((time_vector_4sec(1,i)),1) == k
            x_axis_Labels(i,1) = k;
            k = k + 0.5;
        end
    end

axe2 = subplot(3,1,2);

h3 = heatmap(prepost_dFF_discrete);
%h3.Title = 'dF/F, two seconds pre-/post-bout transition';
h.XLabel = 'Time (seconds)';
h3.YLabel = 'Bout';
h3.Colormap = parula;
h3.GridVisible = 'off';
h3.ColorMethod = 'none';
h3.MissingDataLabel = '';
h3.XDisplayLabels = x_axis_Labels;
h3.YDisplayLabels = y_axis_Labels;

% plot averaged pre/post trace with SEM shading

line3 = zeros(size(prepost_dFF_discrete, 2), 1);

axes3 = subplot(3,1,3);
stdshade(prepost_dFF_discrete, 0.10, 'k', time_vector_4sec);
hold on 
plot(time_vector_4sec, line3, 'k:');
xline(0);
hold off
%title('Averaged dF/F across bouts, first 2 seconds', 'FontSize',10);
ylabel('dF/F');
xlabel('Time (seconds)');

%% Freezing bout-by-bout analysis

%So far, have photometry signal in the array within_Freezing_dFF_discrete which is an m x n matrix where m = the max length of the entire recording session and n = the number of bouts detected. I then
%filled the array such that each row represents the dFF values for a bout and the rest is just NaNs.

%First, let's just output some parameters for each bout in columnar format

summary_Bouts = NaN(num_Bouts, 8);
summary_prepost = NaN(size(prepost_dFF_discrete,1),9);

for i = 1:num_Bouts
    
    %summary parameters for all bouts
    summary_Bouts(i,1) = nanmean(within_Freezing_dFF_discrete(i,:));
    summary_Bouts(i,2) = nanmax(within_Freezing_dFF_discrete(i,:));
    summary_Bouts(i,3) = nanmin(within_Freezing_dFF_discrete(i,:));
    summary_Bouts(i,4) = nanstd(within_Freezing_dFF_discrete(i,:));
    summary_Bouts(i,5) = ((stop_indices(i) - start_indices(i)) / data.streams.g2_2.fs); %the length of the bout in indices/sampling frequency)
    
    % summary parameters for bouts > 2 seconds
    summary_Bouts(i,6) = nanmean(within_Freezing_dFF_discrete_sorted(i,:));
    summary_Bouts(i,7) = nanmax(within_Freezing_dFF_discrete_sorted(i,:));
    summary_Bouts(i,8) = nanmin(within_Freezing_dFF_discrete_sorted(i,:));
    summary_Bouts(i,9) = nanstd(within_Freezing_dFF_discrete_sorted(i,:));    
    summary_Bouts(i,10) = ((stop_indices(i) - start_indices(i)) / data.streams.g2_2.fs);    

end

% Pre/post analyses

for i = 1:size(prepost_dFF_discrete,1)
    
    %parameters for post (0 - 2 seconds)
    summary_prepost(i,1) = nanmean(prepost_dFF_discrete(i, (1:round(2 * data.streams.g2_2.fs))));
    summary_prepost(i,2) = nanmax(prepost_dFF_discrete(i, (1:round(2 * data.streams.g2_2.fs))));
    summary_prepost(i,3) = nanmin(prepost_dFF_discrete(i, (1:round(2 * data.streams.g2_2.fs))));
    summary_prepost(i,4) = nanstd(prepost_dFF_discrete(i, (1:round(2 * data.streams.g2_2.fs))));

    
    %parameters for post (2 - 4 seconds)
    summary_prepost(i,5) = nanmean(prepost_dFF_discrete(i, (round(2 * data.streams.g2_2.fs):round(4 * data.streams.g2_2.fs))));
    summary_prepost(i,6) = nanmax(prepost_dFF_discrete(i, (round(2 * data.streams.g2_2.fs):round(4 * data.streams.g2_2.fs))));
    summary_prepost(i,7) = nanmin(prepost_dFF_discrete(i, (round(2 * data.streams.g2_2.fs):round(4 * data.streams.g2_2.fs))));
    summary_prepost(i,8) = nanstd(prepost_dFF_discrete(i, (round(2 * data.streams.g2_2.fs):round(4 * data.streams.g2_2.fs))));

    
    %timestamp dFF values from beginning, transition, and end of bout
    summary_prepost(i,9) = nanmean(prepost_dFF_discrete(i, 1:round(0.1 * data.streams.g2_2.fs)));
    summary_prepost(i,10) = nanmean(prepost_dFF_discrete(i, (round(1.9 * data.streams.g2_2.fs):round(2.1 * data.streams.g2_2.fs))));
    summary_prepost(i,11) = nanmean(prepost_dFF_discrete(i, (round(3.9 * data.streams.g2_2.fs):round(4 * data.streams.g2_2.fs))));    
    
% below is old code defining start/end as a singular dFF value from the first/last indice of the array.
% instead, I am going to average the first and last 100 ms (shown above), basically equivalent to a 10X downsampling
%     summary_prepost(i,7) = prepost_dFF_discrete(i,1);
%     summary_prepost(i,8) = prepost_dFF_discrete(i,round(size(prepost_dFF_discrete,2)));
%     summary_prepost(i,9) = prepost_dFF_discrete(i,size(prepost_dFF_discrete,2));

end

disp('analysis finished')
%% Export figures

saveas(figure(1), strcat('MATLAB output/Plots/', sheetname, ' - Figure 1.tif'))
saveas(figure(2), strcat('MATLAB output/Plots/', sheetname, ' - Figure 2.tif'))
saveas(figure(3), strcat('MATLAB output/Plots/', sheetname, ' - Figure 3.tif'))
saveas(figure(4), strcat('MATLAB output/Plots/', sheetname, ' - Figure 4.tif'))

end


