% Master organizer file to archive all mice/sessions, calls photometry analysis pipeline as a function
% to iterate through specified combinations of mice and output plots/data in icKO vs. control format

%% To-do list

% Right now, the code successfully iterates through string arrays containing the info for specific sessions. It can output saved figures for each session but
% still need to figure out how to 1) extract the relevant information into either a .mat file or Excel sheet and 2) work with batch data instead of per session
% (like, working with the entire range of dFF values for multiple sessions at once rather than just means, maxs, etc.

% First step is probably to load various single row arrays of dFF values (CFC, HC, freezing, mobile). Maybe have the function return several row of values back to 
% a master matrix here? Probably need a 3D matrix for bout data.

% 3D matrices work but I get 0 values? Figure out why and/or how to get rid of them.
% Preallocate my matrices so it's faster and can specify NaN matrices instead.

%% Housekeeping
    close all; clear all; clc;

%% Databases

% Array structure:
% Rows = one session
% Columns, in order = BLOCKPATH, BLOCKPATH_HC, onset_FF, filename, sheetname
                  
database_cond_ctrl = {'269_simon-181101-144059', '590', 'Photometry #1', '269_cond'; 
    
                      '577_simon-181101-123405', '13', 'Photometry #1', '577_cond'; 

                      '194-181205-171817', '13', 'Photometry #2', '194_cond'; 
            
                      '195-181205-162853', '3', 'Photometry #2', '195_cond'; 
                
                      '197-181205-151205', '12', 'Photometry #2', '197_cond'; 
            
                      '198-181205-153755', '18', 'Photometry #2', '198_cond'; 
                  
                      '271-181205-165220', '10', 'Photometry #2', '271_cond';} ;
                      
                  
database_cond_icKO = {'193-181205-174146', '10', 'Photometry #2', '193_cond'; 
    
                      '199-181205-160529', '19', 'Photometry #2', '199_cond';
                  
                      '501_simon-181101-142602', '4', 'Photometry #1', '504_cond';} ; 
                  

database_24hr_ctrl = {'269_simon-181102-131525', '269_simon-181102-130313', '7', 'Photometry #1', '269_24 hour'; 
    
                      '577_simon-181102-120527', '577_simon-181102-114800', '16', 'Photometry #1', '577_24 hour'; 

                      '194-181206-143038', '194-181206-142217', '20', 'Photometry #2', '194_24 hour'; 
            
                      '195-181206-155238', '195-181206-154500', '6', 'Photometry #2', '195_24 hour'; 
                
                      '197-181206-134625', '197-181206-133923', '154', 'Photometry #2', '197_24 hour'; 
            
                      '198-181206-132458', '198-181206-131704', '11', 'Photometry #2', '198_24 hour'; 
                  
                      '271-181206-140934', '271-181206-140228', '3', 'Photometry #2', '271_24 hour';} ;
                      
                 
database_24hr_icKO = {'193-181206-152955', '193-181206-152111', '8', 'Photometry #2', '193_24 hour'; 
    
                      '199-181206-130406', '199-181206-125646', '6', 'Photometry #2', '199_24 hour';
                  
                      '504_simon-181102-123823', '504_simon-181102-122855', '17', 'Photometry #1', '504_24 hour';} ; 
                  
                  
database_30d_ctrl = {'269_simon-181129-152406', '269_simon-181129-153210', '16', 'Photometry #1', '269_30 Day'; 
    
                      '577_simon-181129-145435', '577_simon-181129-144301', '5', 'Photometry #1', '577_30 day'; 

                      '194-190103-141950', '194-190103-141020', '14', 'Photometry #2', '194_30 day'; 
            
                      '195-190103-135332', '195-190103-134406', '177', 'Photometry #2', '195_30 day'; 
                
                      '197-190103-120954', '197-190103-120125', '14', 'Photometry #2', '197_30 day'; 
            
                      '198-190103-130202', '198-190103-125422', '9', 'Photometry #2', '198_30 day'; 
                  
                      '271-190103-132917', '271-190103-132126', '32', 'Photometry #2', '271_30 day';} ;
                  
                  
database_30d_icKO = {'193-190103-144011', '193-190103-143246', '10', 'Photometry #2', '193_30 day'; 
    
                      '199-190103-123948', '199-190103-123109', '11', 'Photometry #2', '199_30 day';
                            
                      '504_simon-181129-151425', '504_simon-181129-150631', '13', 'Photometry #1', '504_30 day';}  ;                    


%% Iterate through analyses and grab dFF + summary parameters for each session

set(0,'DefaultFigureVisible','off') % suppress plots, doesn't actually work

run_dB = vertcat(database_24hr_ctrl, database_24hr_icKO, database_30d_ctrl, database_30d_icKO); %insert databases here to run, concatenate if necessary
%database_24hr_ctrl, database_24hr_icKO, database_30d_ctrl, database_30d_icKO

sessions = cell(size(run_dB,1),1); 

for i = 1:(size(run_dB,1)) % populate session names
    sessions(i) = run_dB(i,5);
end

run_dB = string(run_dB); %convert character array to cell array of strings

for i = 1:(size(run_dB,1))
    
    BLOCKPATH_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB(i,1))); %TDT import function requires cell or character array input // vhpc_t_PFC_NAc_social_bar-180202-135210
    BLOCKPATH_HC_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB(i,2)));
    onset_FF_input = run_dB(i,3);
    filename_input = run_dB(i,4);
    sheetname_input = run_dB(i,5);

    [m1,m2,m3,m4,m5,m6] = photometry_f(BLOCKPATH_input,BLOCKPATH_HC_input,onset_FF_input,filename_input,sheetname_input); %grab outputs in temporary matrices
    
    %generate 2D or 3D arrays for batch data
    for j = 1:size(m1,2)
        batch_dFF(i,j) = m1(j);
    end
    
    for j = 1:size(m2,2)
        batch_dFF_HC(i,j) = m2(j);
    end
    
    for j = 1:size(m3,2)
        for k = 1:size(m3, 1)
            batch_prepost_dFF(k,j,i) = m3(k,j); %3D array structure: rows = bouts, columns = single dFF values, pages = sessions
        end
    end
 
    for j = 1:size(m4,2)
        batch_summary_Mean(i,j) = m4(j);
    end
    
    for j = 1:size(m5,2)
        for k = 1:size(m5, 1)
            batch_summary_Bouts(k,j,i) = m5(k,j); %3D array structure: rows = bouts, columns = parameters, pages = sessions
        end
    end
    
    for j = 1:size(m6,2)
        for k = 1:size(m6, 1)
            batch_summary_prepost(k,j,i) = m6(k,j); %3D array structure: rows = bouts, columns = parameters, pages = sessions
        end
    end

end   

disp('parse finished')

%% Export Excel sheets

%Export session mean data
t1 = table(sessions, batch_summary_Mean(:,1), batch_summary_Mean(:,2), batch_summary_Mean(:,3), batch_summary_Mean(:,4), batch_summary_Mean(:,5), batch_summary_Mean(:,6), batch_summary_Mean(:,7), batch_summary_Mean(:,8), batch_summary_Mean(:,9), batch_summary_Mean(:,10), batch_summary_Mean(:,11), batch_summary_Mean(:,12), batch_summary_Mean(:,13), batch_summary_Mean(:,14), batch_summary_Mean(:,15), batch_summary_Mean(:,16));
t1.Properties.VariableNames = {'Session', 'CFCmean', 'CFCmin', 'CFCmax', 'CFCSTD', 'HCmean', 'HCmin', 'HCmax', 'HCSTD', 'FreezeMean', 'FreezeMin', 'FreezeMax', 'FreezeSTD', 'MobileMean', 'MobileMin', 'MobileMax', 'MobileSTD'}

writetable(t1, 'MATLAB output/batch_summary_Mean.xls');

%Export sorted bout data

for i = 1:size(batch_summary_Bouts,3)
    
    t2 = table(batch_summary_Bouts(:,1,i), batch_summary_Bouts(:,2,i), batch_summary_Bouts(:,3,i), batch_summary_Bouts(:,4,i), batch_summary_Bouts(:,5,i), batch_summary_Bouts(:,6,i), batch_summary_Bouts(:,7,i), batch_summary_Bouts(:,8,i), batch_summary_Bouts(:,9,i), batch_summary_Bouts(:,10,i));
    t2.Properties.VariableNames = {'Mean', 'Min', 'Max', 'Length', 'STD', 'Mean2s', 'Min2s', 'Max2s', 'Length2s', 'STD2s'}
    writetable(t2, 'MATLAB output/batch_summary_Bouts.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));

end

%Export pre/post analysis

for i = 1:size(batch_summary_prepost,3)
    
    t3 = table(batch_summary_prepost(:,1,i), batch_summary_prepost(:,2,i), batch_summary_prepost(:,3,i), batch_summary_prepost(:,4,i), batch_summary_prepost(:,5,i), batch_summary_prepost(:,6,i), batch_summary_prepost(:,7,i), batch_summary_prepost(:,8,i), batch_summary_prepost(:,9,i), batch_summary_prepost(:,10,i), batch_summary_prepost(:,11,i));
    t3.Properties.VariableNames = {'PreMean', 'PreMax', 'PreMin', 'PreSTD', 'PostMean', 'PostMax', 'PostMin', 'PostSTD', 'StartDFF', 'TransDFF', 'EndDFF'}
    writetable(t3, 'MATLAB output/batch_summary_prepost.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));

end

% %Export peaks analysis
% 
% % summary_peaks just go into a sheet in the overall CFC spreadsheet
%     t4 = table(sessions, batch_summary_Peaks(:,1), batch_summary_Peaks(:,2), batch_summary_Peaks(:,3), batch_summary_Peaks(:,4), batch_summary_Peaks(:,5), batch_summary_Peaks(:,6), batch_summary_Peaks(:,7), batch_summary_Peaks(:,8), batch_summary_Peaks(:,9), batch_summary_Peaks(:,10), batch_summary_Peaks(:,11), batch_summary_Peaks(:,12), batch_summary_Peaks(:,13), batch_summary_Peaks(:,14), batch_summary_Peaks(:,15), batch_summary_Peaks(:,16), batch_summary_Peaks(:,17), batch_summary_Peaks(:,18), batch_summary_Peaks(:,19), batch_summary_Peaks(:,20)); 
%     t4.Properties.VariableNames = {'Session', 'NumPeaksCFC', 'MeanHeightCFC', 'MeanWidthCFC', 'MaxHeightCFC', 'MaxWidthCFC', 'NumPeaksHC', 'MeanHeightHC', 'MeanWidthHC', 'MaxHeightHC', 'MaxWidthHC', 'NumPeaksFRZ', 'MeanHeightFRZ', 'MeanWidthFRZ', 'MaxHeightFRZ', 'MaxWidthFRZ', 'NumPeaksMOB', 'MeanHeightMOB', 'MeanWidthMOB', 'MaxHeightMOB', 'MaxWidthMOB'}
%     writetable(t4, 'MATLAB output/batch_peaks_CFC.xls', 'FileType', 'spreadsheet', 'Sheet', 'summary_Peaks');
% 
% for i = 1:size(batch_peaks_CFC,3)
%     
%     t5 = table(batch_peaks_CFC(:,1,i), batch_peaks_CFC(:,2,i), batch_peaks_CFC(:,3,i), batch_peaks_CFC(:,4,i));
%     t5.Properties.VariableNames = {'Height', 'Loc', 'Width', 'Prominence'}
%     writetable(t5, 'MATLAB output/batch_peaks_CFC.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));
% end
% 
% for i = 1:size(batch_peaks_HC,3)
%     
%     t6 = table(batch_peaks_HC(:,1,i), batch_peaks_HC(:,2,i), batch_peaks_HC(:,3,i), batch_peaks_HC(:,4,i));
%     t6.Properties.VariableNames = {'Height', 'Loc', 'Width', 'Prominence'}
%     writetable(t6, 'MATLAB output/batch_peaks_HC.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));
% end
% 
% for i = 1:size(batch_peaks_Frz,3)
%     
%     t7 = table(batch_peaks_Frz(:,1,i), batch_peaks_Frz(:,2,i), batch_peaks_Frz(:,3,i), batch_peaks_Frz(:,4,i));
%     t7.Properties.VariableNames = {'Height', 'Loc', 'Width', 'Prominence'}
%     writetable(t7, 'MATLAB output/batch_peaks_Frz.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));
% end
% 
% for i = 1:size(batch_peaks_Mob,3)
%     
%     t8 = table(batch_peaks_Mob(:,1,i), batch_peaks_Mob(:,2,i), batch_peaks_Mob(:,3,i), batch_peaks_Mob(:,4,i));
%     t8.Properties.VariableNames = {'Height', 'Loc', 'Width', 'Prominence'}
%     writetable(t8, 'MATLAB output/batch_peaks_Mob.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));
% end

disp('batch finished')

%% Single batch analysis

% Run these sections optionally after the above section

%% Across session averaged CFC trace (with optional HC trace overlay)

line = zeros(1,size(batch_dFF,2));
time_vector_300s = linspace(0, 300, size(batch_dFF, 2));

figure(1);
stdshade(batch_dFF, 0.10, 'k', time_vector_300s);
hold on 
%stdshade(batch_dFF_HC, 0.10, 'k', time_vector_300s); % HC overlay would be nice but I probably messed up the HC for this experiment. I would also have to trim it down to exactly 300s to overlay.
plot(time_vector_300s, line, 'k:');
hold off
%title('Averaged dF/F across bouts, first 2 seconds', 'FontSize',10);
ylabel('dF/F');
xlabel('Time (seconds)');
ylim([-0.2 0.2]);

%% Across session averaged pre/post trace

% allocate space for (x * z) * y array
plot_prepost = NaN(((size(batch_prepost_dFF,1)) * (size(batch_prepost_dFF,3))), size(batch_prepost_dFF,2));

for i = 1:size(batch_prepost_dFF,3)    
    plot_prepost = vertcat(plot_prepost, batch_prepost_dFF(:,:,i));
end

plot_prepost(any(isnan(plot_prepost),2)|any(plot_prepost==0,2),:) = []; % trim zero and NaN-only rows
line = zeros(1,size(plot_prepost,2));
time_vector_4sec = linspace(-2,2,size(plot_prepost,2));

figure(1);
stdshade(plot_prepost, 0.10, 'k', time_vector_4sec);
hold on 
plot(time_vector_4sec, line, 'k:');
xline(0);
hold off
%title('Averaged dF/F across bouts, first 2 seconds', 'FontSize',10);
ylabel('dF/F');
xlabel('Time (seconds)')
%ylim([-0.06 0.06]);

%% Two batch analysis

% Run this section to load data into two arrays to make comparisons between groups

set(0,'DefaultFigureVisible','off') % suppress plots, doesn't actually work

run_dB = vertcat(database_30d_ctrl); %insert databases here to run, concatenate if necessary
run_dB2 = vertcat(database_30d_icKO); %insert another set of sessions for comparison
%database_24hr_ctrl, database_24hr_icKO, database_30d_ctrl, database_30d_icKO

sessions = cell(size(run_dB,1),1); 
sessions2 = cell(size(run_dB2,1),1); 

for i = 1:(size(run_dB,1)) % populate session names
    sessions(i) = run_dB(i,5);
end

for i = 1:(size(run_dB2,1)) % populate session names
    sessions2(i) = run_dB2(i,5);
end

run_dB = string(run_dB); %convert character array to cell array of strings
run_dB2 = string(run_dB2); 

for i = 1:(size(run_dB,1))
    
    BLOCKPATH_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB(i,1))); %TDT import function requires cell or character array input // vhpc_t_PFC_NAc_social_bar-180202-135210
    BLOCKPATH_HC_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB(i,2)));
    onset_FF_input = run_dB(i,3);
    filename_input = run_dB(i,4);
    sheetname_input = run_dB(i,5);

    [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11] = photometry_f(BLOCKPATH_input,BLOCKPATH_HC_input,onset_FF_input,filename_input,sheetname_input); %grab outputs in temporary matrices
    
    %generate 2D or 3D arrays for batch data
    for j = 1:size(m1,2)
        batch_dFF(i,j) = m1(j);
    end
    
    for j = 1:size(m2,2)
        batch_dFF_HC(i,j) = m2(j);
    end
    
    for j = 1:size(m3,2)
        for k = 1:size(m3, 1)
            batch_prepost_dFF(k,j,i) = m3(k,j); %3D array structure: rows = bouts, columns = single dFF values, pages = sessions
        end
    end
 
    for j = 1:size(m4,2)
        batch_summary_Mean(i,j) = m4(j);
    end
    
    for j = 1:size(m5,2)
        for k = 1:size(m5, 1)
            batch_summary_Bouts(k,j,i) = m5(k,j); %3D array structure: rows = bouts, columns = parameters, pages = sessions
        end
    end
    
    for j = 1:size(m6,2)
        for k = 1:size(m6, 1)
            batch_summary_prepost(k,j,i) = m6(k,j); %3D array structure: rows = bouts, columns = parameters, pages = sessions
        end
    end
    
    % import batch peaks data
    for j = 1:size(m7, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m7,1)
            batch_peaks_CFC(j,k,i) = m7(k,j);
        end
    end
    
    for j = 1:size(m8, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m8,1)
            batch_peaks_HC(j,k,i) = m8(k,j);
        end
    end
    
    for j = 1:size(m9, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m9,1)
            batch_peaks_Frz(j,k,i) = m9(k,j);
        end
    end

    for j = 1:size(m10, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m10,1)
            batch_peaks_Mob(j,k,i) = m10(k,j);
        end
    end
    
    for j = 1:size(m11,2);
        batch_summary_Peaks(i,j) = m11(j);
    end
       
end   

%run second parse

for i = 1:(size(run_dB2,1))
    
    BLOCKPATH_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB2(i,1))); %TDT import function requires cell or character array input // vhpc_t_PFC_NAc_social_bar-180202-135210
    BLOCKPATH_HC_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB2(i,2)));
    onset_FF_input = run_dB2(i,3);
    filename_input = run_dB2(i,4);
    sheetname_input = run_dB2(i,5);

    [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11] = photometry_f(BLOCKPATH_input,BLOCKPATH_HC_input,onset_FF_input,filename_input,sheetname_input); %grab outputs in temporary matrices
    
    %generate 2D or 3D arrays for batch data
    for j = 1:size(m1,2)
        batch_dFF2(i,j) = m1(j);
    end
    
    for j = 1:size(m2,2)
        batch_dFF_HC2(i,j) = m2(j);
    end
    
    for j = 1:size(m3,2)
        for k = 1:size(m3, 1)
            batch_prepost_dFF2(k,j,i) = m3(k,j); %3D array structure: rows = bouts, columns = single dFF values, pages = sessions
        end
    end
 
    for j = 1:size(m4,2)
        batch_summary_Mean2(i,j) = m4(j);
    end
    
    for j = 1:size(m5,2)
        for k = 1:size(m5, 1)
            batch_summary_Bouts2(k,j,i) = m5(k,j); %3D array structure: rows = bouts, columns = parameters, pages = sessions
        end
    end
    
    for j = 1:size(m6,2)
        for k = 1:size(m6, 1)
            batch_summary_prepost2(k,j,i) = m6(k,j); %3D array structure: rows = bouts, columns = parameters, pages = sessions
        end
    end
    
    % import batch peaks data
    for j = 1:size(m7, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m7,1)
            batch_peaks_CFC2(j,k,i) = m7(k,j);
        end
    end
    
    for j = 1:size(m8, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m8,1)
            batch_peaks_HC2(j,k,i) = m8(k,j);
        end
    end
    
    for j = 1:size(m9, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m9,1)
            batch_peaks_Frz2(j,k,i) = m9(k,j);
        end
    end

    for j = 1:size(m10, 2)  %3D array structure: rows = peaks, columns = parameters, pages = sessions
        for k = 1:size(m10,1)
            batch_peaks_Mob2(j,k,i) = m10(k,j);
        end
    end
    
    for j = 1:size(m11,2);
        batch_summary_Peaks2(i,j) = m11(j);
    end
       
end   

% now, have across session raw dFF values in various format and summary data for overall session/bouts/peaks
% can do statistics within Matlab to compare two groups, overlay group averaged traces, etc


%% Conditioning analysis

% If running conditioning sessions, run this section separately and not the above two sections

set(0,'DefaultFigureVisible','off') % suppress plots, doesn't actually work

run_dB = vertcat(database_cond_icKO, database_cond_ctrl); %insert databases here to run, concatenate if necessary
sessions = cell(size(run_dB,1),1); % store session names for tabular output

for i = 1:(size(run_dB,1)) % populate session names
    sessions(i) = run_dB(i,4);
end

run_dB = string(run_dB); %convert character array to cell array of strings

for i = 1:(size(run_dB,1))
    
    BLOCKPATH_input = cellstr(strcat('vhpc_t_PFC_NAc_social_bar-180202-135210/', run_dB(i,1))); %TDT import function requires cell or character array input // vhpc_t_PFC_NAc_social_bar-180202-135210
    onset_FF_input = run_dB(i,2);
    filename_input = run_dB(i,3);
    sheetname_input = run_dB(i,4);
    
    [m1,m2,m3] = photometry_f_cond(BLOCKPATH_input,onset_FF_input, sheetname_input, filename_input); %grab outputs in temporary matrices
    
    % populate aligned pre/post dFFs across sessions
    for j = 1:size(m1,2)
        for k = 1:3
            batch_cond_dFF(k,j,i) = m1(k,j); %3D array structure: rows = shocks (3), columns = dFF, pages = sessions
        end
    end
    
    % populate mean across shock parameters
    for j = 1:size(m2,2)
        batch_summary_cond(i,j) = m2(j);
    end
    
    % populate peri-shock parameters
    for j = 1:size(m3,2)
        for k = 1:3
            batch_summary_prepostShock(k,j,i) = m3(k,j); %3D array structure: rows = shocks (3), columns = parameters, pages = sessions
        end
    end
end

%% Batch analysis (conditioning)

% Plot shock response for first shock as heatmap, and averaged shock responses
% code here is kinda just freeform, not well commented

% First, convert 3D matrix to 2D matrix, convert 3D pages = 2D rows

batch_cond_dFF_shock1 = NaN(size(batch_cond_dFF,3),size(batch_cond_dFF,2));
%%

for i = 1:size(batch_cond_dFF,3)
    batch_cond_dFF_shock1(i,:) = batch_cond_dFF(1,:,i);
end

%%
clf
figure(1);

time_vector_5to5 = linspace(-5, 5, size(batch_cond_dFF,2)); 
x_axis_Labels = NaN(size(time_vector_5to5,2), 1);
y_axis_Labels = NaN(5,1);

k = -5;
    for i = 1:size(time_vector_5to5,2)
        if round((time_vector_5to5(1,i)),1) == k
            x_axis_Labels(i,1) = k;
            k = k + 5;
        end
    end

h = heatmap(squeeze(batch_cond_dFF_shock1));
%h.Title = '5 seconds pre-/post-shock dF/F';
h.XLabel = 'Time (seconds)';
h.YLabel = '';
h.Colormap = parula;
h.GridVisible = 'off';
h.ColorMethod = 'none';
h.MissingDataLabel = '';
h.XDisplayLabels = x_axis_Labels;
h.YDisplayLabels = y_axis_Labels;
h.FontName = 'Helvetica'
h.FontSize = 20

%%

figure(2);

axes = gca;
batch_cond_dFF_shocks = NaN((2*size(batch_cond_dFF,3)),size(batch_cond_dFF,2));

r = 1;

% populate shock1
  for i = 1:size(batch_cond_dFF,3)
    batch_cond_dFF_shocks(r,:) = batch_cond_dFF(1,:,i);
    r = r + 1;
  end

% populate shock3
  for i = 1:size(batch_cond_dFF,3)
    batch_cond_dFF_shocks(r,:) = batch_cond_dFF(3,:,i);
    r = r + 1;
  end

line2 = zeros(size(batch_cond_dFF_shocks, 2), 1);

stdshade(batch_cond_dFF_shocks, 0.10, 'k', time_vector_5to5);
hold on 
%v = vline(0, 'k');
plot(time_vector_5to5, line2, 'k:');
hold off
%title('Averaged pre-/post-shock dF/F', 'FontSize',10);
ylabel('');
xlabel('Time (seconds)');
axes.FontSize = 20;


%% Export Excel sheets (conditioning analyses)

%Export across-shock mean/peak
t1 = table(sessions, batch_summary_cond(:,1), batch_summary_cond(:,2), batch_summary_cond(:,3), batch_summary_cond(:,4), batch_summary_cond(:,5), batch_summary_cond(:,6), batch_summary_cond(:,7), batch_summary_cond(:,8), batch_summary_cond(:,9), batch_summary_cond(:,10), batch_summary_cond(:,11), batch_summary_cond(:,12));
t1.Properties.VariableNames = {'Session', 'PreshockMean', 'PostShock1Mean', 'PostShock2Mean', 'PostShock3Mean', 'PreshockMax', 'PostShock1Max', 'PostShock2Max', 'PostShock3Max', 'PreshockSTD', 'PostShock1STD', 'PostShock2STD', 'PostShock3STD'}

writetable(t1, 'MATLAB output/batch_conditioning.xls', 'FileType', 'spreadsheet', 'Sheet', 'summary_cond');

%Export sorted bout data

for i = 1:size(batch_summary_prepostShock,3)
    
    t2 = table(batch_summary_prepostShock(:,1,i), batch_summary_prepostShock(:,2,i), batch_summary_prepostShock(:,3,i), batch_summary_prepostShock(:,4,i), batch_summary_prepostShock(:,5,i));
    t2.Properties.VariableNames = {'PreMean', 'PostMean', 'PreMax', 'PostMax', 'Peak'}
    writetable (t2, 'MATLAB output/batch_conditioning.xls', 'FileType', 'spreadsheet', 'Sheet', char(cellstr(sessions(i,1))));

end

disp('batch finished')
