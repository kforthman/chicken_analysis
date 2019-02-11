


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

trial = {'T1','T2'};

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

k=1;
for i = 1:length(subjids)
    
    disp(i)
    
    if availableData{i,2} && availableData{i,3}
        for j = 1:2
            filename = [pathFrom 'wave2-CT-' subjids{i} '-' trial{j} '.csv'];
            
            fid = fopen(filename);
            header = strsplit(fgetl(fid), ',');
            %colNames = strsplit(fgetl(fid), ',');
            fclose(fid);
            
            ua_pos = find(strcmp(header, 'UserAGENT'));
            useragent = strjoin(header((ua_pos+1):end), ',');
            header = header(1:(ua_pos+1));
            header((ua_pos+1)) = cellstr(useragent);
            
            header = reshape(header, 2, []);
            header = cell2struct(header(2,:), header(1,:), 2);
            
            data = readtable(filename, 'ReadVariableNames', true, 'Delimiter', ',', 'HeaderLines', 1);
            
            % ID
            splt = strsplit(header.Orginal_File_Name, {'-','.'});
            id{k} = splt{2};
            
            % Pattern
            pattern{k} = header.Pattern;
            
            t{k} = trial{j};
            
            k = k + 1;
        end
    end
end

id = id.';
t = t.';
pattern = pattern.';

myTable = table(id,t,pattern);
writetable(myTable, fullfile(pathTo, 'blockPatterns.csv'))