function [datastruct] = createStruct(filename)
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
id = splt{2};

% Trial
trial = splt{3};

% Version:
version = header.Version;

% Type:
type = header.Type;

% Pattern:
pattern = header.Pattern;

if strcmp(pattern, '1')
    block_order = [3, 1, 5, 6, 2, 4];
elseif strcmp(pattern, '2')
    block_order = [6, 2, 4, 3, 1, 5];
elseif strcmp(pattern, '3')
    block_order = [1, 3, 5, 2, 4, 6];
elseif strcmp(pattern, '4')
    block_order = [2, 4, 6, 1, 3, 5];
end

% N: total trials for the session
practice_location = strcmp(data.trial_type, 'practice');
main_location = strcmp(data.trial_type, 'main');

N = sum(main_location);
data = data(main_location,:);

% signaled: whether the true mean was revealed or not at the end of each trial
if strcmp(type, 'estimate')
    signaled = ones(N,1);
elseif strcmp(type, 'predict')
    signaled = zeros(N,1);
end

% sigma
% based on block type
sigma = zeros(N,1);
for j = 1:N
    eval(['sigma(' num2str(j) ') = str2num(header.block' num2str(data.block_num(j)) '_sigma);'])
end


% muall: raw [x,y] coordinates of the two chickens (left and right). Center of the screen = left chicken X + ((right chicken X - left chicken X)/2)
muall = [str2num(header.chicken_left_x),str2num(header.chicken_left_y);...
    str2num(header.chicken_right_x),str2num(header.chicken_right_y)];

% midpt: The midpoint between the two chickens.
midpt = mean(muall(:,1));

% X: [x,y] coordinates of each star on each trial
X = [data.egg_x_position, data.egg_y_position];


% muinds: The actual triangle that generated the star (1 for left, 2 for right)
muinds = zeros(N,1);
for j = 1:N
    if [data.response(j), data.result(j)] == [1,1]
        muinds(j) = 1;
    elseif [data.response(j), data.result(j)] == [2,1]
        muinds(j) = 2;
    elseif [data.response(j), data.result(j)] == [1,2]
        muinds(j) = 2;
    elseif [data.response(j), data.result(j)] == [2,2]
        muinds(j) = 1;
    end
end

% musgn: Signaled chicken from previous trial,
% (neg => left, pos => right):
musgn = sign(muall(muinds,1)-midpt);
% If the correct chicken was not revealed at the end of the trial, musgn is set to 0.
musgn(signaled==0) = 0;
musgn = [0;musgn(1:end-1)];

% pred: the triangle the star appears to have come from (not which one it actually came from)
pred = data.response;

% rt: reaction times
rt = data.response_time_sec;

% H: current hazard rate
% based on block type
H = zeros(N,1);
for j = 1:N
    eval(['H(' num2str(j) ') = str2num(header.block' num2str(data.block_num(j)) '_hazard);'])
end

%H = H.';

% Hset: gives the hazard rate of each block.
Hset = [H(100) H(300) H(500) H(700) H(900) H(1100)];

% Sset: gives the sigma of each block
Sset = [sigma(100) sigma(300) sigma(500) sigma(700) sigma(900) sigma(1100)];

% cp: changepoint. 1 if the mean is different from the previous trial. 0 otherwise.
cp = zeros(N,1);
cp(1) = 1;
for j = 2:N
    if muinds(j) ~= muinds(j-1)
        cp(j) = 1;
    end
end



% r: trial since last change-point
r = zeros(N,1);
r(1) = 1;
for j = 2:N
    if cp(j) == 1
        r(j) = 1;
    else
        r(j) = r(j-1) + 1;
    end
end

% Hset: Set of hazard rates for the current session (either two hazard rates - 1000 trials each, or 5 hazard rates - 400 trials each)
% Appears to not be needed

% cp: whether or not the star appears to have come from the same triangle as the previous trial. Note that this is not the ACTUAL change-point of the triangle
% Appears to not be needed

% r: trial since last change-point as defined above
% Appears to not be needed

%tau
% Appears to not be needed
datastruct = cell(6,1);
block_num = table2array(data(:,'block_num'));
for j = 1:6
session = block_num == block_order(j);

datastruct{j} = struct('ID', id, 'trial', trial, 'block', block_order(j), 'N', sum(session), 'sigma', Sset(j), 'muall', muall, 'X', X(session, :), ...
    'pred', pred(session), 'muinds', muinds(session), 'rt', rt(session), 'H', H(session), ...
    'signaled', signaled(session), 'cp', cp(session), 'r', r(session),'Hset', Hset(j));
end

