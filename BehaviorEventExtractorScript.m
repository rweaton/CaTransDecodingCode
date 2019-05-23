% Set path to file that contains data
PathToFile = 'G:\SwapFolder_retracking\20181205\20181205_both_001_1-03Dat.csv';

% Read in table from csv file
B = readtable(PathToFile);

% Identify "pre-pause" delay in timestamp record
Filt = [false; diff(B.Timestamp) > (B.Timestamp(end) - B.Timestamp(end - 1))];
Indices = transpose(1:length(Filt));
StartIndex = Indices(Filt);

% Subtract out pre-pause delay on startup
B.Timestamp = B.Timestamp(:) - B.Timestamp(StartIndex);

% Remove rows from table corresponding to pause time.
B = B(B.Timestamp >= 0.,:);

% Generate a list of indices for the resulting filtered table
Indices = transpose(1:length(B.Timestamp));

% All right-hand target 0 enter events
rh_t0_enters = B.Timestamp([false; diff(B.EV1_1) == 1]);

% All right-hand center exit events
rh_ct_exits = B.Timestamp([false; diff(B.EV1_5) == -1]);

% All left-hand target 1 enter events
lh_t1_enters = B.Timestamp([false; diff(B.EV1_9) == 1]);

% All left-hand center exit events
lh_ct_exits = B.Timestamp([false; diff(B.EV1_11) == -1]);

% Select nearest center target exit that precedes each target entry
rh_ct_exits_filtered = NaN*ones(1, length(rh_t0_enters));
lh_ct_exits_filtered = NaN*ones(1, length(lh_t1_enters));

% Set the maxiumum time allowed for a center target exit to precede a
% peripheral target entrance.
move_time_limit = 1.00;  % units of seconds.

% Iterate through peripheral target enter events to find nearest preceding
% center target exit event.
rh_ct_exits_indices = 1:length(rh_ct_exits);

for i = 1:length(rh_t0_enters)
    
    % Select preceding exits only!
    PrecedeFilt = (rh_ct_exits < rh_t0_enters(i));
    
    % Acquire subset of indices that point timestamps that precede current
    % event timestamp.
    PrecedeIndices = rh_ct_exits_indices(PrecedeFilt);
    
    % Find nearest timestamp
    [dist, index] = min(rh_t0_enters(i) - rh_ct_exits(PrecedeIndices));
    
    % Include only if it lies within the limit of movement time
    if dist <= move_time_limit
        
        rh_ct_exits_filtered(i) = rh_ct_exits(PrecedeIndices(index));
    end
    
end

% Now repeat for left hand
lh_ct_exits_indices = 1:length(lh_ct_exits);

for i = 1:length(lh_t1_enters)
    
    % Select preceding exits only!
    PrecedeFilt = (lh_ct_exits < lh_t1_enters(i));
    
    % Acquire subset of indices that point timestamps that precede current
    % event timestamp.
    PrecedeIndices = lh_ct_exits_indices(PrecedeFilt);
    
    % Find nearest timestamp
    [dist, index] = min(lh_t1_enters(i) - lh_ct_exits(PrecedeIndices));
    
    % Include only if it lies within the limit of movement time
    if dist <= move_time_limit
        
        lh_ct_exits_filtered(i) = lh_ct_exits(PrecedeIndices(index));
    end
    
end

all_events_new = {};
all_events_new {1,1} = 'rh_ct_exits_filtered'; all_events_new{1,2} = transpose(rh_ct_exits_filtered);
all_events_new{2,1} = 'lh_ct_exits_filtered'; all_events_new{2,2} = transpose(lh_ct_exits_filtered);



