function fitAdaptivityModel(filename)
% function fitAdaptivityModel
%
% Fit testing data to the Adaptivity Model
%
%
% Altered so that pct (percent correct) and H estimates are computed independently for each
% value of sigma. There is no parameter for learning rate.
%
% Edited by Katherine Forthman 08.20.2018

%% Initialize values, parameters are:
%   1. H_subjective
%   2. noise in the decision variable (DV)
%   3. lapse rate
% These are the parameters that will be estimated using the optimization
% function later in the program.
%-kf-% J denotes the log prior odds. All instances of J have been converted to H.
MU_DIST = 150; %-kf-% Distance between each source (in pixels)
SIGMAS  = round([110 160]); %-kf-% The standard deviation of the star distribution around triangle (also pixels)
nsigmas = length(SIGMAS);  %-kf-% The number of different sigma values
ub      = [1 10 .5];       %-kf-% upper boundaries on parameter estimates
lb      = [0 10e-8 10e-8]; %-kf-% lower boundaries on parameter estimates
inits   = [0.5 .2 .0025];  %-kf-% initial parameter estimates

%% Get information about data files
%-kf-% The function getDataInfo() is called to download the following
%-kf-% objects:
%-kf-% file_list := list of filenames of each participant's data
%-kf-% analysis_data_dir := The directory of the Analysis folder
%-kf-% raw_data_dir := The directory of the Raw folder
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;
subID = arrayfun(@(x) ...
    subsref(...
    strsplit(file_list{x},{'.','_'}),...
    struct('type','()','subs',{{1}})...
    ), 1:length(file_list)).';

%% Loop through the subjects
num_subjects = length(file_list);
%-kf-% The following empty objects are created to be filled in the following loop.
%-kf-% fits := contains the parameter estimates for each subject
%-kf-% LLRs := The log likelihood ratio provided by the star position
%               (sensory evidence) for each participant.
%-kf-% pcts := percent correct
fits         = nan(num_subjects, size(inits, 2));
LLRs         = nan(num_subjects, 1);
pcts         = nan(num_subjects, 1); % percent correct

for ss = 1:num_subjects
    
    disp(['Simulation ' num2str(ss)])
    
    % get the data file
    data_filename = fullfile(raw_data_dir, file_list{ss});
    load(data_filename)
    
    %-kf-% x is star position
    % collect some data for each session
    % in the structure session_data:
    %   1. likelihood of x given left
    %   2. likelihood of x given right
    %   3. choice
    %   4. H %-kf-% the log of the posterior odds
    %   5. Signaled mean from previous trial (if given)
    ctmp = nan(1, 3); % sigma, # correct, N
    
    
    data = datastruct;
    
    
    %-kf-% Center of screen (x coordinate):
    midpt = mean(data.muall(:,1));
    %-kf-% Signaled chicken from previous trial,
    %-kf-% (neg => left, pos => right):
    musgn = sign(data.muall(data.muinds,1)-midpt);
    musgn = [0;musgn(1:end-1)];
    
    % collect data in a single matrix
    % session_data is cell array (per session), each cell has a single
    %   matrix with rows as trials, columns are:
    %   1. likelihood of x given left
    %   2. likelihood of x given right
    %   3. choice
    %   4. H
    %   5. Signaled mean from previous trial (if given)
    session_data = cat(2, ...
        normpdf(data.X(:,1)-midpt, MU_DIST/2, data.sigma), ...
        normpdf(data.X(:,1)-midpt,-MU_DIST/2, data.sigma), ...
        double(data.pred==2));%, ... %-kf-% data.pred is the participant choice (1 for left, 2 for right) -> (0 for left, 1 for right)
    %data.H(:,1), ...
    %musgn);
    
    % save stuff to compute % correct per sigma
    ctmp = [ ...
        data.sigma, ...
        sum(data.pred==data.muinds), ...
        size(data.pred,1)];
    
    % do the fit for this subject
    myFun = @(x)fitAdaptivityModel_err(x, session_data);
    
    % now... fit it
    options = optimoptions(@fmincon,'Algorithm','interior-point');
    fits(ss,:) = ...
        fmincon(myFun, inits, [], [], [], [], lb, ub, [], options);
    
    
    N{ss} = data.N;
    H_true{ss} = data.Htrue;
    H_subj{ss} = data.Hsubj;
    sigma{ss}  = data.sigma;
    pct{ss}    = (ctmp(2)/ctmp(3))*100;
end



H_subj_EST = fits(:,1);
noise_in_DV = fits(:,2);
lapse_rate = fits(:,3);

% Save it
if nargin<1 || (nargin>=1 && ~ischar(filename))
    filename1 = 'adaptivityModelFits.csv';
    filename2 = 'adaptivityModelFits.mat';
    filename3 = 'percentCorrect.csv';
end

N = N.';
H_subj = H_subj.';
H_true = H_true.';
sigma = sigma.';
pct = pct.';

finalTable = table(subID,N,sigma,H_true,H_subj,H_subj_EST,noise_in_DV,lapse_rate,pct);
writetable(finalTable, fullfile(analysis_data_dir, filename1))
