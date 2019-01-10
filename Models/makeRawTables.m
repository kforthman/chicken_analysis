


% Designate the destination for the formated data
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;
num_sessions = 4;
num_subjects = length(file_list);
subID = arrayfun(@(x) subsref(strsplit(file_list{x},{'.','_'}),struct('type','()','subs',{{2}})), 1:length(file_list)).';

for ss = 1:num_subjects
    
    disp(ss)
    
    % get the data file
    data_filename = fullfile(raw_data_dir, file_list{ss});
    load(data_filename)
    
    for ii = 1:num_sessions
        
        % get the appropriate data
        eval(['data=dataT' num2str(ii) ';'])
        names = {'egg_x';'egg_y';'correct_choice';'H_true';'sigma';'prediction';'reaction_time'};
        m = length(data.H);
        finalTable = table(data.X(:,1),...
            data.X(:,2),...
            data.muinds,...
            data.H,...
            ones(m,1)*data.sigma,...
            data.pred,...
            data.rt, ...
            'VariableNames',names);
        
        filename1 = [subID{ss} '_T' num2str(ii) '.csv'];
        writetable(finalTable, [raw_data_dir '/' filename1])
        
    end
end