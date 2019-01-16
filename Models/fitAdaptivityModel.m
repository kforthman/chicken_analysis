function fitAdaptivityModel(filename)
% function fitAdaptivityModel
%
% Fit testing data to the Adaptivity Model
%
% This takes <1 min on a MacBook Pro 2.8 GHz Intel Core i7
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line
%     learning in an unpredictable environment"
%
% Altered so that pct (percent correct) is computed independently for each
% value of sigma.
%
% Edited by Katherine Forthman 08.20.2018

%% Initialize values, params are:
%-kf-% J denotes the log prior odds
%   1. J_default
%   2. m_H
%   3. noise in the DV %-kf-% discrete values?
%   4. lapse rate
MU_DIST = 150; %-kf-% Distance between each source
SIGMAS  = round([33 140]); %-kf-% The standard deviation of the star distribution around triangle
nsigmas = length(SIGMAS); %-kf-%
ub      = [ 10  10 10 .02]; %-kf-% upper boundaries on parameter estimates
lb      = [-10 -10  0 .0001]; %-kf-% lower boundaries on parameter estimates
inits   = [0.45 0.4 0.3 .0025]; %-kf-%  initial parameter estimates ???

%% Get information about data files
%-kf-% The function getDataInfo() is called to download the following
%-kf-% objects:
%-kf-% file_list := list of filenames of each participant's data
%-kf-% analysis_data_dir := The directory of the Analysis folder
%-kf-% raw_data_dir := The directory of the Raw folder
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;

%% Loop through the subjects
num_subjects = length(file_list);
%-kf-% The following empty objects are created for future use.
%-kf-% fits := contains the parameter estimates for each subject?
%-kf-% LLRs := The log likelihood ratio provided by the star position
%(sensory evidence) for each participant.
%-kf-% pcts := accuracy per sigma
fits         = nan(num_subjects, size(inits, 2));
LLRs         = nan(num_subjects, 1);
pcts_33      = nan(num_subjects, 1); % percent correct per sigma
pcts_140     = nan(num_subjects, 1);

for ss = 1:num_subjects
    
    disp(ss)
    
    % get the data file
    data_filename = fullfile(raw_data_dir, file_list{ss});
    load(data_filename)
    
    %-kf-% x is star position
    % collect some data for each session:
    %   1. likelihood of x given left
    %   2. likelihood of x given right
    %   3. choice
    %   4. logJ0 = log[H/(1-H)] %-kf-% the log of the posterior odds?
    %   5. Signaled mean from previous trial (if given)
    num_sessions = 4;
    session_data = cell(num_sessions, 1);
    ctmp = nan(num_sessions, 3); % sigma, # correct, N
    
    for ii = 1:num_sessions
        
        % get the appropriate data
        eval(['data=dataT' num2str(ii) ';'])
        
        % useful stuff
        %-kf-% Center of screen (x coordinate):
        midpt = mean(data.muall(:,1));
        %-kf-% Signaled chicken from previous trial,
        %-kf-% (neg => left, pos => right):
        musgn = sign(data.muall(data.muinds,1)-midpt);
        %-kf-% If the true mean was not revealed at the end of the trial, musgn is set to 0.
        musgn(data.signaled==0) = 0;
        
        musgn = [0;musgn(1:end-1)];
        
        % collect data in a single matrix
        % session_data is cell array (per session), each cell has a single
        %   matrix with rows as trials, columns are:
        %   1. likelihood of x given left
        %   2. likelihood of x given right
        %   3. choice
        %   4. logJ0 = log[H/(1-H)]
        %   5. Signaled mean from previous trial (if given)
        session_data{ii} = cat(2, ...
            normpdf(data.X(:,1)-midpt, MU_DIST/2, data.sigma), ...
            normpdf(data.X(:,1)-midpt,-MU_DIST/2, data.sigma), ...
            double(data.pred==2), ... %-kf-% data.pred is the actual triangle that generated the star (1 for left, 2 for right)
            log(data.H(:,1)./(1-data.H(:,1))), ... % data.H(:,1), ...
            musgn);
        
        % save stuff to compute % correct per sigma
        ctmp(ii,:) = [ ...
            data.sigma, ...
            sum(data.pred==data.muinds), ...
            size(data.pred,1)];
    end
    
    % do the fit for this subject
    myFun = @(x)fitAdaptivityModel_err(x, session_data);
    
    % now... fit it
    [fits(ss,:), LLRs(ss)] = ...
        fmincon(myFun, inits, [], [], [], [], lb, ub, []);
    % optimoptions(@fmincon,'Algorithm','interior-point'));
    
    % compute overall accuracy per sigma
    Lsig = round(ctmp(:,1)) == SIGMAS(1);
    if any(Lsig)
        pcts_33(ss) = sum(ctmp(Lsig,2))./sum(ctmp(Lsig,3)).*100;
    end
    Lsig = round(ctmp(:,1)) == SIGMAS(2);
    if any(Lsig)
        pcts_140(ss) = sum(ctmp(Lsig,2))./sum(ctmp(Lsig,3)).*100;
    end
end



J_default = fits(:,1);
m_H = fits(:,2);
noise_in_DV = fits(:,3);
lapse_rate = fits(:,4);

H_default = arrayfun(@(x) 1/(1 + exp(-J_default(x))), 1:length(J_default)).';

subID = arrayfun(@(x) subsref(strsplit(file_list{x},{'.','_'}),struct('type','()','subs',{{2}})), 1:length(file_list)).';

% Save it
if nargin<1 || (nargin>=1 && ~ischar(filename))
    filename1 = 'adaptivityModelFits.csv';
    filename2 = 'adaptivityModelFits.mat';
end

finalTable = table(subID,H_default,m_H,noise_in_DV,lapse_rate,LLRs,pcts_33,pcts_140);
writetable(finalTable, fullfile(analysis_data_dir, filename1))
save(fullfile(analysis_data_dir, filename2),'LLRs','pcts_33','pcts_140','H_default','m_H','noise_in_DV','lapse_rate','subID','J_default');
