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
%   1. H_subjective
%   3. noise in the decision variable (DV)
%   4. lapse rate
MU_DIST = 150; %-kf-% Distance between each source
SIGMAS  = round([33 140]); %-kf-% The standard deviation of the star distribution around triangle
nsigmas = length(SIGMAS); %-kf-%
ub      = [1 10 .5]; %-kf-% upper boundaries on parameter estimates
lb      = [0 10e-8 10e-8]; %-kf-% lower boundaries on parameter estimates
inits   = [0.5 .2 .0025]; %-kf-%  initial parameter estimates

%% Get information about data files
%-kf-% The function getDataInfo() is called to download the following
%-kf-% objects:
%-kf-% file_list := list of filenames of each participant's data
%-kf-% analysis_data_dir := The directory of the Analysis folder
%-kf-% raw_data_dir := The directory of the Raw folder
[file_list, analysis_data_dir, raw_data_dir] = getDataInfo;
subID = arrayfun(@(x) subsref(strsplit(file_list{x},{'.','_'}),struct('type','()','subs',{{2}})), 1:length(file_list)).';

%% Loop through the subjects
num_subjects = length(file_list);
%-kf-% The following empty objects are created for future use.
%-kf-% fits := contains the parameter estimates for each subject?
%-kf-% LLRs := The log likelihood ratio provided by the star position
%(sensory evidence) for each participant.
%-kf-% pcts := accuracy per sigma
num_sessions = 12;
fits         = nan(num_subjects*num_sessions, size(inits, 2));
LLRs         = nan(num_subjects*num_sessions, 1);
pcts_33      = nan(num_subjects, 1); % percent correct per sigma
pcts_140     = nan(num_subjects, 1);

for kk = 0:0.1:1
    eval(['fits_' num2str(kk*10) '= nan(num_subjects*num_sessions, size(inits, 2))'])
end

for ss = 1:num_subjects
    
    disp(['Participant ' num2str(ss)])
    
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
    
    session_data = cell(num_sessions, 1);
    ctmp = nan(num_sessions, 3); % sigma, # correct, N
    
    for ii = 1:num_sessions
        
        % get the appropriate data
        if ii<=6
            eval(['data=dataT1B' num2str(ii) ';'])
            disp(['Participant ' num2str(ss) ', Session T1B' num2str(ii)])
        else
            eval(['data=dataT2B' num2str(ii-6) ';'])
            disp(['Participant ' num2str(ss) ', Session T2B' num2str(ii-6)])
        end
        
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
        %   4. H
        %   5. Signaled mean from previous trial (if given)
        session_data{ii} = cat(2, ...
            normpdf(data.X(:,1)-midpt, MU_DIST/2, data.sigma), ...
            normpdf(data.X(:,1)-midpt,-MU_DIST/2, data.sigma), ...
            double(data.pred==2));%, ... %-kf-% data.pred is the participant choice (1 for left, 2 for right) -> (0 for left, 1 for right)
        %data.H(:,1), ...
        %musgn);
        
        % save stuff to compute % correct per sigma
        ctmp(ii,:) = [ ...
            data.sigma, ...
            sum(data.pred==data.muinds), ...
            size(data.pred,1)];
        
        % do the fit for this subject
        myFun = @(x)fitAdaptivityModel_err(x, session_data{ii});
        
        % now... fit it
        this.index = (ss-1)*12+ii;
        options = optimoptions(@fmincon,'Algorithm','interior-point');
        for kk = 0:0.1:1
            disp(kk)
            inits   = [kk .2 .0025];
            eval(['fits_' num2str(kk*10) '(this.index,:) = fmincon(myFun, inits, [], [], [], [], lb, ub, [], options);'])
        end
        
        ID{this.index}     = data.ID;
        trial{this.index}  = data.trial;
        block{this.index}  = num2str(data.block);
        H_true{this.index} = data.Hset;
        sigma{this.index}  = data.sigma;
        pct{this.index}    = (ctmp(ii,2)/ctmp(ii,3))*100;
    end
    
    
    
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


for kk = 0:0.1:1
    ind = num2str(kk*10);
    eval(['H_subjective_' ind ' = fits_' ind '(:,1);'])
    eval(['noise_in_DV_' ind ' = fits_' ind '(:,2);'])
    eval(['lapse_rate_' ind ' = fits_' ind '(:,3);'])
end


% Save it
if nargin<1 || (nargin>=1 && ~ischar(filename))
    filename1 = 'adaptivityModelFits_stabilityTest.csv';
    filename2 = 'adaptivityModelFits_stabilityTest.mat';
    filename3 = 'percentCorrect_stabilityTest.csv';
end

ID = ID.';
trial = trial.';
block = block.';
H_true = H_true.';
sigma = sigma.';
pct = pct.';

finalTable = table(ID,trial,block,sigma,H_true,...
    H_subjective_0,H_subjective_1,H_subjective_2,H_subjective_3,H_subjective_4,H_subjective_5,H_subjective_6,H_subjective_7,H_subjective_8,H_subjective_9,H_subjective_10,...
    noise_in_DV_0, noise_in_DV_1 ,noise_in_DV_2 ,noise_in_DV_3 ,noise_in_DV_4 ,noise_in_DV_5 ,noise_in_DV_6 ,noise_in_DV_7 ,noise_in_DV_8 ,noise_in_DV_9 ,noise_in_DV_10 ,...
    lapse_rate_0,  lapse_rate_1,  lapse_rate_2,  lapse_rate_3,  lapse_rate_4,  lapse_rate_5,  lapse_rate_6,  lapse_rate_7,  lapse_rate_8,  lapse_rate_9,  lapse_rate_10,   ...
    pct);
correctTable = table(subID,pcts_33,pcts_140);
%writetable(finalTable, fullfile(analysis_data_dir, filename1))
%writetable(correctTable, fullfile(analysis_data_dir, filename3))
%save(fullfile(analysis_data_dir, filename2),'ID','trial','block','H_subjective','noise_in_DV','lapse_rate');
