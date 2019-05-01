function [datastruct] = createStruct_Sim(N, sigma, Htrue, Hsubj, index)
% N: total trials for the session
% Htrue: current hazard rate
% sigma

%---------------------------------------------------------------------------------------------------------
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
    if muinds(i) == 1
        x(i) = random('Normal', muall(2,1), sigma);
        y(i) = random('Normal', muall(2,2), sigma);
    else
        x(i) = random('Normal', muall(1,1), sigma);
        y(i) = random('Normal', muall(1,2), sigma);
    end
end

% Now, based on the data above, the model will make a prediction of which chicken the egg came from.
l = zeros(N,1);
psi = zeros(N,1);
llr = zeros(N,1);
ll = @(this_x, this_mu, this_sigma) 0.5*(this_x - this_mu)^2/this_sigma^2 + 0.5*log(this_sigma^2);
l(1) = 0;
psi(1) = 0;
llr(1) = 0;
for n = 2:N
    psi(n) = l(n - 1) + ...
        log((1 - Hsubj)/Hsubj + exp(-l(n - 1))) - ...
        log((1 - Hsubj)/Hsubj + exp( l(n - 1)));
    
    llr(n) = -2*...
        (...
        ll(x(n), muall(2,1), sigma) -...
        ll(x(n), muall(1,1), sigma)...
        );
    
    l(n) = psi(n) + llr(n);
end
% pred: the participant's response.
pred = (l < 0)+1; % for estimation task
pred2 = (psi < 0)+1; % for prediction task
crct = (muinds == pred); % for estimation task
crct2 = (muinds == pred2); % for prediction task

% X: [x,y] coordinates of each egg on each trial
X = [x, y];

% Appears to not be needed
datastruct = struct('N', N, 'sigma', sigma, 'muall', muall, 'X', X, ...
    'pred', pred2, 'muinds', muinds, 'Htrue', Htrue, 'Hsubj', Hsubj, 'cp', swch);

%save(['Data/Simulated_Raw/sim_' num2str(N) '_' num2str(sigma) '_' num2str(Htrue) '_' num2str(Hsubj) '_' num2str(index) '.mat'],'datastruct');

save(sprintf('Data/Simulated_Raw/sim_%d_%d_%.2f_%.2f_%d.mat', N, sigma, Htrue, Hsubj, index),'datastruct')
end

