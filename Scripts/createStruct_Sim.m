function [datastruct] = createStruct_Sim(sigma, Htrue, Hsubj, lapse, noise, index)
% N: total trials for the session
% Htrue: current hazard rate
% sigma

%---------------------------------------------------------------------------------------------------------
N = 150;
% muall: raw [x,y] coordinates of the two chickens (left and right). Center of the screen = left chicken X + ((right chicken X - left chicken X)/2)
muall = [-75,0;75,0];
left_coord = [-75,0];
right_coord = [75,0];

left_val = 1;
right_val = 2;

% Then, generate some random data to simulate the behavior of the eggs in each trial.
r = random('Uniform',0,1,[N,1]); % create a string of random values between 0 and 1
swch = (r <= Htrue); % select those values greater than hTrue to get switch points.
% muinds: The actual chicken that generated the egg (1 for left, 2 for right)
muinds = zeros(N,1); % initialize muinds, the vector that indicates the correct chicken
val = 1;
for i = 1:N
    if swch(i) == 1 % If switch
        val = (val == 1)+1; % value changes
    end
    muinds(i) = val; % muinds stores the correct answer
end

% Generate the x,y position based on the correct chicken
x = zeros(N,1);
y = zeros(N,1);
for i = 1:N
    if muinds(i) == right_val
        x(i) = random('Normal', right_coord(1), sigma);
        y(i) = random('Normal', right_coord(2), sigma);
    else
        x(i) = random('Normal', left_coord(1), sigma);
        y(i) = random('Normal', left_coord(2), sigma);
    end
end

%%
%   1. likelihood of x given right
l_right = normpdf(x, right_coord(1), sigma);
%   2. likelihood of x given left
l_left = normpdf(x, left_coord(1), sigma);
% LLRn represents the sensory evidence. LLRn is the log likelihood ratio
% provided ny the star position on that trial.
llr = log(l_right./l_left);

%%
% Now, based on the data above, the model will make a prediction of which chicken the egg came from.
l = zeros(N,1);
psi = zeros(N,1);

l(1) = 0;
psi(1) = 0;

for n = 2:N
    psi(n) = l(n - 1) + ...
        log((1 - Hsubj)/Hsubj + exp(-l(n - 1))) - ...
        log((1 - Hsubj)/Hsubj + exp( l(n - 1)));
    
    l(n) = psi(n) + llr(n);
end
% pred: the participant's response.
% pred = (l > 0)+1; % for estimation task
% pred2 = (psi > 0)+1; % for prediction task
% crct = (muinds == pred); % for estimation task
% crct2 = (muinds == pred2); % for prediction task

%predhat: the participant's response with lapse and noise.
choicehat = lapse+(1-2*lapse)./(1+exp(-l/noise)); % for estimation task
choice2hat = lapse+(1-2*lapse)./(1+exp(-psi/noise)); % for prediction task

pred  = binornd(1,  choicehat)+1; % for estimation task
pred2 = binornd(1, choice2hat)+1; % for prediction task

crct = (muinds == pred); % for estimation task
crct2 = (muinds == pred2); % for prediction task

% X: [x,y] coordinates of each egg on each trial
X = [x, y];


%save(['Data/Simulated_Raw/sim_' num2str(N) '_' num2str(sigma) '_' num2str(Htrue) '_' num2str(Hsubj) '_' num2str(index) '.mat'],'datastruct');

str_sigma = num2str(sigma, '%03.0f');
str_hTrue = num2str(Htrue*100, '%02.0f');
str_hSubj = num2str(Hsubj*100, '%02.0f');
str_lapse = num2str(lapse*100, '%02.0f');
str_noise = num2str(noise*100, '%03.0f');
str_index = num2str(index, '%02.0f');

datastruct = struct('N', 50, 'sigma', sigma, 'noise', noise, 'lapse', lapse, 'muall', muall, 'X', X(1:50,:), ...
    'pred', pred2(1:50), 'muinds', muinds(1:50), 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch(1:50), 'correct', crct2(1:50));
save(...
    ['../Data/Simulated_nl_Raw/sim-'...
    'N_050' ...
    '-sigma_' str_sigma ...
    '-hT_'    str_hTrue ...
    '-hS_'    str_hSubj ...
    '-lapse_' str_lapse ...
    '-noise_' str_noise ...
    '-' str_index ...
    '.mat'],'datastruct')

datastruct = struct('N', 100, 'sigma', sigma, 'noise', noise, 'lapse', lapse, 'muall', muall, 'X', X(1:100,:), ...
    'pred', pred2(1:100), 'muinds', muinds(1:100), 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch(1:100), 'correct', crct2(1:100));
save(...
    ['../Data/Simulated_nl_Raw/sim-'...
    'N_100' ...
    '-sigma_' str_sigma ...
    '-hT_'    str_hTrue ...
    '-hS_'    str_hSubj ...
    '-lapse_' str_lapse ...
    '-noise_' str_noise ...
    '-' str_index ...
    '.mat'],'datastruct')

datastruct = struct('N', 150, 'sigma', sigma, 'noise', noise, 'lapse', lapse, 'muall', muall, 'X', X, ...
    'pred', pred2, 'muinds', muinds, 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch, 'correct', crct2);
save(...
    ['../Data/Simulated_nl_Raw/sim-'...
    'N_150' ...
    '-sigma_' str_sigma ...
    '-hT_'    str_hTrue ...
    '-hS_'    str_hSubj ...
    '-lapse_' str_lapse ...
    '-noise_' str_noise ...
    '-' str_index ...
    '.mat'],'datastruct')
end

