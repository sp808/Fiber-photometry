%% Function declaration

% analysis for conditioning sessions only

function [prepost_shock_dFF, summary_cond, summary_prepostShock] = photometry_f_cond(BLOCKPATH_input,onset_FF_input, sheetname_input, filename_input)

%% Manually import session information

% Use this section if not running in batch function mode, comment out otherwise
close all; clear all; clc;

BLOCKPATH = 'vhpc_t_PFC_NAc_social_bar-180202-135210/195-181205-162853'
onset_FF = 3
sheetname = '195_cond'
filename = 'Photometry #2'

%% Import data and define manual inputs
% 
% BLOCKPATH = char(BLOCKPATH_input)
% onset_FF = str2double(onset_FF_input) %Manually input
% sheetname = char(sheetname_input)
% filename = char(filename_input) % Photometry #2 - 193_21 day - 271_30 day - Event Logs // Photometry #1 - 269_24 hour - 577_30 day - Event Logs

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

offset_FF = onset_FF + 420; %start time is manually input, stop time is 5 min later

% [min_value_onset, start_idx] = min(abs(time_g22 - onset_FF)); %subtract element-wise by the timestamped value and then find the index of the lowest value (closest to zero)
% [min_value_offset, stop_idx] = min(abs(time_g22 - offset_FF));

start_idx = ceil(onset_FF * data.streams.g2_2.fs);
stop_idx = ceil(offset_FF * data.streams.g2_2.fs);

dFF_trim = NaN(1, round(420 * data.streams.g2_2.fs));

k = 1;
for i = start_idx:stop_idx
    dFF_trim(k) = dFF(i); % align trimmed trace to t = 0
    k = k + 1;
end

%% Parse baseline and post-shock intervals 

% convert shock times to indices
shock1 = round(300 * data.streams.g2_2.fs);
shock2 = round(330 * data.streams.g2_2.fs);
shock3 = round(360 * data.streams.g2_2.fs);

summary_cond = NaN(1,8);

summary_cond(1) = nanmean(dFF(1:shock1));
summary_cond(2) = nanmean(dFF(shock1:shock2));
summary_cond(3) = nanmean(dFF(shock2:shock3));
summary_cond(4) = nanmean(dFF(shock3:(round(offset_FF * data.streams.g2_2.fs))));

summary_cond(5) = nanmax(dFF(1:shock1));
summary_cond(6) = nanmax(dFF(shock1:shock2));
summary_cond(7) = nanmax(dFF(shock2:shock3));
summary_cond(8) = nanmax(dFF(shock3:(round(offset_FF * data.streams.g2_2.fs))));

summary_cond(9) = nanstd(dFF(1:shock1));
summary_cond(10) = nanstd(dFF(shock1:shock2));
summary_cond(11) = nanstd(dFF(shock2:shock3));
summary_cond(12) = nanstd(dFF(shock3:(round(offset_FF * data.streams.g2_2.fs))));

%% Prepost analysis

% 5 seconds pre, 5 seconds post, pre/post means and peak response
% so, first construct a 3 x 2 array to hold the start/stop indices for each of the three 10 sec peri-shock intervals
% then, populate a 3 row array of length (10 sec * sampling frequency) to hold the dFF values for each peri-shock interval
% extract parameters and plots

%populate array with indices for 10 second peri-shock intervals
shock_prepost_indices = [round((onset_FF + 295) * data.streams.g2_2.fs), round((onset_FF + 305) * data.streams.g2_2.fs);
                         round((onset_FF + 325) * data.streams.g2_2.fs), round((onset_FF + 335) * data.streams.g2_2.fs);
                         round((onset_FF + 355) * data.streams.g2_2.fs), round((onset_FF + 365) * data.streams.g2_2.fs);];

%populate array with indices for 2 second actual shock intervals (for plot visualization)
shock_visual_indices = [round((onset_FF + 300) * data.streams.g2_2.fs), round((onset_FF + 302) * data.streams.g2_2.fs);
                        round((onset_FF + 330) * data.streams.g2_2.fs), round((onset_FF + 332) * data.streams.g2_2.fs);
                        round((onset_FF + 360) * data.streams.g2_2.fs), round((onset_FF + 362) * data.streams.g2_2.fs);];                     

prepost_shock_dFF = NaN(3, round(10 * data.streams.g2_2.fs));
shock_visual = NaN(1, size(dFF_trim,2)); %visualize NaN vs. 0 on aligned trace plot
summary_prepostShock = NaN(3, 5);

% populate prepost_shock_dFF for each of the three shocks and align t = 0
for i = 1:3
    k = 1;
    for j = shock_prepost_indices(i,1) : shock_prepost_indices(i,2)
        prepost_shock_dFF(i,k) = dFF(j);
        k = k +1;
    end
    
    %also populate shock_visual for plotting while we're here
    for j = shock_visual_indices(i,1) : shock_visual_indices(i,2)
        shock_visual(1,(j - floor(onset_FF * data.streams.g2_2.fs))) = 0;
    end
end

% extract pre mean, post mean, pre max, pre mean, and overall peak
for i = 1:3
    summary_prepostShock(i,1) = nanmean(prepost_shock_dFF(i, 1:round(5 * data.streams.g2_2.fs)));
    summary_prepostShock(i,2) = nanmean(prepost_shock_dFF(i, round(5 * data.streams.g2_2.fs):round(10 * data.streams.g2_2.fs)));
    summary_prepostShock(i,3) = nanmax(prepost_shock_dFF(i, 1:round(5 * data.streams.g2_2.fs)));
    summary_prepostShock(i,4) = nanmax(prepost_shock_dFF(i, round(5 * data.streams.g2_2.fs):round(10 * data.streams.g2_2.fs)));
    summary_prepostShock(i,5) = nanmax(prepost_shock_dFF(i,:));
end

%% Plot aligned trace with shock annotations and peri-shock heatmaps + averaged traces

figure(1);

time_vector_420sec = linspace(0, 420, size(dFF_trim,2));

ax1 = subplot (3,1,1);
line1 = zeros(1, size(time_vector_420sec,2));
plot(time_vector_420sec, dFF_trim(1,:), 'k');
hold on
plot(time_vector_420sec, line1,'k:');
plot(time_vector_420sec, shock_visual(1,:), 'y', 'LineWidth', 30, 'Color', [0.93,0.69,0.13]);
hold off
%title('dFF over conditioning session', 'FontSize',10);
%ylabel('');
xlim([0 420]);
xlabel('Time (seconds)');
ax1.FontSize = 20;

%% Peri event heatmap

ax2 = subplot (3,1,2);

time_vector_10sec = linspace(0, 10, size(prepost_shock_dFF,2)); 
time_vector_5to5 = linspace(-5, 5, size(prepost_shock_dFF,2)); 
x_axis_Labels = NaN(size(time_vector_10sec,2), 1);

k = -5;
    for i = 1:size(time_vector_10sec,2)
        if round((time_vector_10sec(1,i)),1) == k
            x_axis_Labels(i,1) = k;
            k = k + 5;
        end
    end

h = heatmap(prepost_shock_dFF);
%h.Title = '5 seconds pre-/post-shock dF/F';
%h.XLabel = 'Time (seconds)';
h.YLabel = 'Shock';
h.Colormap = parula;
h.GridVisible = 'off';
h.ColorMethod = 'none';
h.MissingDataLabel = '';
h.XDisplayLabels = x_axis_Labels;

ax3 = subplot (3,1,3);

line2 = zeros(size(prepost_shock_dFF, 2), 1);

stdshade(prepost_shock_dFF, 0.10, 'k', time_vector_5to5);
hold on 
v = vline(0, 'k');
plot(time_vector_5to5, line2, 'k:');
hold off
%title('Averaged pre-/post-shock dF/F', 'FontSize',10);
ylabel('dF/F (%)');
xlabel('Time (seconds)');

%% Export plots

saveas(figure(1), strcat('MATLAB output/Plots/Conditioning', sheetname, ' - Figure 5.tif'))

end
