% Format data from test version 2.2.1 of chicken task to be read into
% fitAdaptivityModel.m (developed by Glaze et al.)
%
% Created by Katherine Forthman, 09.17.2018

% Designate the file in which the data is located
pathFrom = '/Volumes/T1000/BehavioralTasks/MTURK/ChickenTask/data/tasks/data_for_katie/';

% Designate the destination for the formated data
pathTo = '/Users/kclary/Dropbox (LIBR)/02_Projects/Chicken Task/Chicken_code copy/Data/';

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
        
%         dataT1B1 = struct('N', 199, 'sigma', dataT1.sigma, 'muall', dataT1.muall, 'X', dataT1.X(1:199, :), ...
%             'pred', dataT1.pred(1:199), 'muinds', dataT1.muinds(1:199), 'rt', dataT1.rt(1:199), 'H', dataT1.H(1:199), ...
%             'signaled', dataT1.signaled(1:199), 'cp', dataT1.cp(1:199), 'r', dataT1.r(1:199),'Hset', dataT1.Hset(1));
%         
%         dataT1B2 = struct('N', 200, 'sigma', dataT1.sigma, 'muall', dataT1.muall, 'X', dataT1.X(200:399, :), ...
%             'pred', dataT1.pred(200:399), 'muinds', dataT1.muinds(200:399), 'rt', dataT1.rt(200:399), 'H', dataT1.H(200:399), ...
%             'signaled', dataT1.signaled(200:399), 'cp', dataT1.cp(200:399), 'r', dataT1.r(200:399),'Hset', dataT1.Hset(2));
%         
%         dataT1B3 = struct('N', 200, 'sigma', dataT1.sigma, 'muall', dataT1.muall, 'X', dataT1.X(400:599, :), ...
%             'pred', dataT1.pred(400:599), 'muinds', dataT1.muinds(400:599), 'rt', dataT1.rt(400:599), 'H', dataT1.H(400:599), ...
%             'signaled', dataT1.signaled(400:599), 'cp', dataT1.cp(400:599), 'r', dataT1.r(400:599),'Hset', dataT1.Hset(3));
%         
%         
%         dataT2B1 = struct('N', 200, 'sigma', dataT2.sigma, 'muall', dataT2.muall, 'X', dataT2.X(1:200, :), ...
%             'pred', dataT2.pred(1:200), 'muinds', dataT2.muinds(1:200), 'rt', dataT2.rt(1:200), 'H', dataT2.H(1:200), ...
%             'signaled', dataT2.signaled(1:200), 'cp', dataT2.cp(1:200), 'r', dataT2.r(1:200),'Hset', dataT2.Hset(1));
%         
%         dataT2B2 = struct('N', 200, 'sigma', dataT2.sigma, 'muall', dataT2.muall, 'X', dataT2.X(201:400, :), ...
%             'pred', dataT2.pred(201:400), 'muinds', dataT2.muinds(201:400), 'rt', dataT2.rt(201:400), 'H', dataT2.H(201:400), ...
%             'signaled', dataT2.signaled(201:400), 'cp', dataT2.cp(201:400), 'r', dataT2.r(201:400),'Hset', dataT2.Hset(2));
%         
%         dataT2B3 = struct('N', 200, 'sigma', dataT2.sigma, 'muall', dataT2.muall, 'X', dataT2.X(401:600, :), ...
%             'pred', dataT2.pred(401:600), 'muinds', dataT2.muinds(401:600), 'rt', dataT2.rt(401:600), 'H', dataT2.H(401:600), ...
%             'signaled', dataT2.signaled(401:600), 'cp', dataT2.cp(401:600), 'r', dataT2.r(401:600),'Hset', dataT2.Hset(3));
%         
        
        save([pathTo '/Raw/data_' subjids{i} '.mat'], 'dataT1', 'dataT2', 'dataT3', 'dataT4')
    else
        continue
    end
    
end