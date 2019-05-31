function err_ = fitAdaptivityModel_err(params, session_data)
% function err_ = fitAdaptivityModel_err(params, session_data)
%
% params are:
%   1. H_subjective
%   2. noise in the decision variable (DV)
%   3. lapse rate
%
% session_data is cell array (per session), each cell has a single
%   matrix with rows as trials, columns are:
%   1. likelihood of x given left
%   2. likelihood of x given right
%   3. choice
%   4. H
%   5. Signaled mean from previous trial (if given)
%
% OUTPUT
% data: structure containing model simulations:
%   data.resp: model predicted responses
%   data.Lprior: model prior belief strength
%   data.Pprior: model probability of choosing option 1 (isntead of 0)
%   data.Lpost: model posterior belief strength (after hearing the tone)
%   data.x: observed stimuli
%   data.params: parameters used for simulation
%
%
%
% Modified by Katherine Forthman, 02/20/2019

n = length(session_data(:,1));

%Initialize the data I want to output
resp_sim   = zeros(n,1);
Lprior = NaN(n,1);
Pprior = NaN(n,1);
Lpost  = NaN(n,1);

%Get choice variability and lapses
noise = params(2);
lapse = params(3);

% LLRn represents the sensory evidence. LLRn is the log likelihood ratio
% provided ny the star position on that trial.
LLR = log(session_data(:,1)./session_data(:,2));

% H is altered to simplify the equation.
% it represents (1-H)/H
%H = (1-session_data(:,4))./session_data(:,4);
H = (1-params(1))/params(1);

   %Run simulation
for i=1:n
    if i == 1
        Lprior(i) = 0;
    else
        Lprior(i) = Lpost(i-1)+log(H+exp(-Lpost(i-1)))-log(H+exp(Lpost(i-1)));
    end
    Lpost(i) = Lprior(i)+LLR(i);
    rprob = lapse+(1-2*lapse)/(1+exp(-Lprior(i)/noise)); %Probability of making a rightward response (1) - computed using prior odds since this is a prediction task

    Pprior(i) = rprob;

    %Simulate response using response probabilities
    if rand < rprob
        resp_sim(i) = 1;
    end
end

%resp = session_data(:,3);
%e = -sum(((1-resp).*log(1-Lpost))+(resp.*log(Lpost)));


   % LogLR from predicted and actual choices
   choicehat = lapse+(1-2*lapse)./(1+exp(-Lprior/noise));
   % choicehat = 1./(1+exp(-Lpost));
   choice    = session_data(:,3);
   e = -sum(((1-choice).*log(1-choicehat))+(choice.*log(choicehat)));
   %e   = -sum((1-choice).*log(1-choicehat)) - sum(choice.*log(choicehat));
   % e = -sum(((1-resp).*log(1-pc))+(resp.*log(pc))); %Cross entropy loss function (same as negative log likelihood)
   % e   = -sum( ((choice == choicehat)-.5)*2 .* log(abs(Lpost)) );

err_ = e;
