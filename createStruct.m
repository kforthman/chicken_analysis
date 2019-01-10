function [datastruct_1, datastruct_2] = createStruct(filename)
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

% Version:
version = header.Version;

% Type:
type = header.Type;

% Pattern:
pattern = header.Pattern;

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

session_1 = (sigma == sigma(1));
session_2 = ~session_1;

% Hset: Set of hazard rates for the current session (either two hazard rates - 1000 trials each, or 5 hazard rates - 400 trials each)
% Appears to not be needed

% cp: whether or not the star appears to have come from the same triangle as the previous trial. Note that this is not the ACTUAL change-point of the triangle
% Appears to not be needed

% r: trial since last change-point as defined above
% Appears to not be needed

%tau
% Appears to not be needed

datastruct_1 = struct('N', sum(session_1), 'sigma', sigma(1), 'muall', muall, 'X', X(session_1, :), ...
    'pred', pred(session_1), 'muinds', muinds(session_1), 'rt', rt(session_1), 'H', H(session_1), ...
    'signaled', signaled(session_1), 'cp', cp(session_1), 'r', r(session_1),'Hset', Hset(1:3));
datastruct_2 = struct('N', sum(session_2), 'sigma', sigma(sum(session_1)+1), 'muall', muall, 'X', X(session_2, :), ...
    'pred', pred(session_2), 'muinds', muinds(session_2), 'rt', rt(session_2), 'H', H(session_2),...
    'signaled', signaled(session_2), 'cp', cp(session_2), 'r', r(session_2),'Hset', Hset(4:6));
end