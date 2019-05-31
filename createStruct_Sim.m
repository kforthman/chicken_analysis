function [datastruct] = createStruct_Sim(sigma, Htrue, Hsubj, index)
% N: total trials for the session
% Htrue: current hazard rate
% sigma

%---------------------------------------------------------------------------------------------------------
N = 150;
% muall: raw [x,y] coordinates of the two chickens (left and right). Center of the screen = left chicken X + ((right chicken X - left chicken X)/2)
muall = [-75,0;75,0];

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
    if muinds(i) == 2
        x(i) = random('Normal', muall(2,1), sigma);
        y(i) = random('Normal', muall(2,2), sigma);
    else
        x(i) = random('Normal', muall(1,1), sigma);
        y(i) = random('Normal', muall(1,2), sigma);
    end
end

%%
%   1. likelihood of x given left
l_left = normpdf(x, muall(2,1)/2, sigma);
%   2. likelihood of x given right
l_right = normpdf(x, muall(1,1)/2, sigma);
% LLRn represents the sensory evidence. LLRn is the log likelihood ratio
% provided ny the star position on that trial.
llr = log(l_left./l_right);

%%
% Now, based on the data above, the model will make a prediction of which chicken the egg came from.
l = zeros(N,1);
psi = zeros(N,1);
%llr = zeros(N,1);
ll = @(this_x, this_mu, this_sigma) 0.5*(this_x - this_mu)^2/this_sigma^2 + 0.5*log(this_sigma^2);
l(1) = 0;
psi(1) = 0;
%llr(1) = 0;
for n = 2:N
    psi(n) = l(n - 1) + ...
        log((1 - Hsubj)/Hsubj + exp(-l(n - 1))) - ...
        log((1 - Hsubj)/Hsubj + exp( l(n - 1)));
    
%     llr(n) = -2*...
%         (...
%         ll(x(n), muall(2,1), sigma) -...
%         ll(x(n), muall(1,1), sigma)...
%         );
    
    l(n) = psi(n) + llr(n);
end
% pred: the participant's response.
pred = (l > 0)+1; % for estimation task
pred2 = (psi > 0)+1; % for prediction task
crct = (muinds == pred); % for estimation task
crct2 = (muinds == pred2); % for prediction task

% X: [x,y] coordinates of each egg on each trial
X = [x, y];


%save(['Data/Simulated_Raw/sim_' num2str(N) '_' num2str(sigma) '_' num2str(Htrue) '_' num2str(Hsubj) '_' num2str(index) '.mat'],'datastruct');

str_m = num2str(N, '%03.0f');
str_sigma = num2str(sigma, '%03.0f');
str_hTrue = num2str(Htrue*100, '%02.0f');
str_hSubj = num2str(Hsubj*100, '%02.0f');
str_index = num2str(index, '%02.0f');

datastruct = struct('N', 50, 'sigma', sigma, 'muall', muall, 'X', X(1:50,:), ...
    'pred', pred2(1:50), 'muinds', muinds(1:50), 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch(1:50));
save(['Data/Simulated_Raw/sim050' str_sigma str_hTrue str_hSubj str_index '.mat'],'datastruct')

datastruct = struct('N', 100, 'sigma', sigma, 'muall', muall, 'X', X(1:100,:), ...
    'pred', pred2(1:100), 'muinds', muinds(1:100), 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch(1:100));
save(['Data/Simulated_Raw/sim100' str_sigma str_hTrue str_hSubj str_index '.mat'],'datastruct')

datastruct = struct('N', 150, 'sigma', sigma, 'muall', muall, 'X', X, ...
    'pred', pred2, 'muinds', muinds, 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch);
save(['Data/Simulated_Raw/sim150' str_sigma str_hTrue str_hSubj str_index '.mat'],'datastruct')
end

