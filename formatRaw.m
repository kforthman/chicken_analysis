% Format data from test version 2.2.1 of chicken task to be read into
% fitAdaptivityModel.m (developed by Glaze et al.)
%
% Created by Katherine Forthman, 09.17.2018

% Designate the file in which the data is located
pathFrom = '/Volumes/T1000/BehavioralTasks/MTURK/ChickenTask/data/tasks/data_for_katie/';

% Designate the destination for the formated data
pathTo = '/Users/kclary/Dropbox (LIBR)/02_Projects/Chicken Task/Chicken_code/Data/';

% Pull the names of all csv files from the raw data directory
rawFilenames = struct2table(dir([pathFrom '*.csv']));

% Split the filenames to separate subject name from session number and
% other naming conventions.
file_details = {};
j=0;
for i = 1:height(rawFilenames)
    nSplit = [rawFilenames.name{i}, strsplit(rawFilenames.name{i},{'-', '.'})];
    if strcmp(nSplit(2),'wave2') && length(nSplit)==6
        file_details(i-j,:) = nSplit;
    else
        j = j+1;
    end
end

% List all the subjects with data in the raw data directory
subjids = unique(file_details(:,4));


% Do the subjects have data for session 1?
T1 = cell(length(subjids),1);
for i = 1:length(subjids)
    T1{i,1} = ismember(array2table([subjids(i), 'T1']), array2table(file_details(:,4:5)), 'rows');
end

% Do the subjects have data for session 2?
T2 = cell(length(subjids),1);
for i = 1:length(subjids)
    T2{i,1} = ismember(array2table([subjids(i), 'T2']), array2table(file_details(:,4:5)), 'rows');
end

% A table of each subject and their sessions completed.
availableData = cell2table([subjids, T1, T2]);

save([pathTo '/Analysis/subj_info.mat'], 'subjids')

% For each subject, create the data structure. The data structure from each
% session is combined into a single file for the subject and saved in the
% output directory. If both sessions are not available, no file is created.
for i = 1:length(subjids)
    
    disp(i)
    
    if availableData{i,2} && availableData{i,3}
        filename = [pathFrom 'wave2-CT-' subjids{i} '-T1.csv'];
        [dataT1, dataT2] = createStruct(filename);

        filename = [pathFrom 'wave2-CT-' subjids{i} '-T2.csv'];
        [dataT3, dataT4] = createStruct(filename);
        
        save([pathTo '/Raw/data_' subjids{i} '.mat'], 'dataT1', 'dataT2', 'dataT3', 'dataT4')
    else
        continue
    end

end