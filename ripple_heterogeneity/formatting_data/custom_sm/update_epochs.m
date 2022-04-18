% populate session epochs by the below logic
%
% a.	Cage1 (animal inside its own home cage, usually on top of the arena to see video) ~1h
% b.	Maze1 (empty arena) ~5min
% c.	Cage2 (animal in its own home cage) ~1h
% d.	Maze2 (arena with empty pencil cups at two opposite corners) ~5min
% e.	Cage3 (animal in its own home cage) ~1h
% f.	Maze3 (arena with two stimulus novel animals, S1 and S2, under each cup) ~5min
% g.	Maze4 (arena with S2 and S1, in opposite positions compare to f) ~5min
% h.	Cage4 (animal in its own home cage) ~1h
% i.	Maze5 (arena with one of the previous stimuli, S1 or S2, and a novel one, N) ~5min

labels = table();
labels.label = {'a','b','c','d','e','f','g','h','i'}';
labels.behavioralParadigm = {'animal in its own home cage',...
    'empty arena',...
    'animal in its own home cage',...
    'arena with empty pencil cups at two opposite corners',...
    'animal in its own home cage',...
    'arena with two stimulus novel animals, S1 and S2, under each cup',...
    'arena with S2 and S1, in opposite positions',...
    'animal in its own home cage',...
    'arena with one of the previous stimuli, S1 or S2, and a novel one, N'}';
labels.environment = {'sleep','box','sleep','box','sleep','box','box','sleep','box'}';
 
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\mouse_sessions.csv');
for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    % check if only single epoch exist and run session template
    if length(session.epochs) == 1
        session = sessionTemplate(basepath,'basename',basename);
        save(fullfile(basepath,[basename,'.session.mat']),'session')
    end
    % check if only single epoch still exist and get older files to help
    if length(session.epochs) == 1
        try
            session = figure_out_epoch_by_behavior_times(basepath,session);
        catch
            % use n samples in dat files
            session = figure_out_epoch_by_dat_files(basepath,basename,session);
        end
    end
    for ep = 1:length(session.epochs)
        if isfield(session.epochs{ep},'environment')
            continue
        end
        label = extractBetween(session.epochs{ep}.name,basename,'_');
        idx = contains(labels.label, label);
        if ~any(idx)
            session.epochs{ep}.environment = 'unknown';
            continue
        end
        session.epochs{ep}.behavioralParadigm = labels.behavioralParadigm{idx};
        session.epochs{ep}.environment = labels.environment{idx};
    end
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end

function session = figure_out_epoch_by_behavior_times(basepath,session)
behavior_times_file = fullfile(basepath,'oldfiles','behavior_times.mat');
if exist(behavior_times_file,'file')
    load(behavior_times_file)
    behavior_times_fields = fields(behavior_times);
    % find subfolders with 'day'
    files = dir(fullfile(basepath,'*day*'));
    sessions = {files([files.isdir]).name}';
    for ep = 1:length(sessions)
        session.epochs{ep}.name = sessions{ep};
        session.epochs{ep}.startTime = behavior_times.(behavior_times_fields{ep})(1);
        session.epochs{ep}.stopTime = behavior_times.(behavior_times_fields{ep})(2);
    end
end
save(fullfile(basepath,[basename,'.session.mat']),'session')
end

function session = figure_out_epoch_by_dat_files(basepath,basename,session)
% find subfolders with 'day'
files = dir(fullfile(basepath,'*day*'));
sessions = {files([files.isdir]).name}';
% start = 0;
for ep = 1:length(sessions)
    fileName = fullfile(basepath,sessions{ep},'amplifier.dat');
    % use dir to get size of dat file in bytes
    filenamestruct = dir(fileName);
    % determine number of bytes per sample
    dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); 
    % Number of samples per channel
    ep_times(ep) = filenamestruct.bytes / (session.extracellular.nChannels*dataTypeNBytes) / session.extracellular.sr;  
end
ep_times = [0,ep_times];
ep_times = cumsum(ep_times);

for ep = 1:length(sessions)-1
    session.epochs{ep}.name = sessions{ep};
    session.epochs{ep}.startTime = ep_times(ep);
    session.epochs{ep}.stopTime = ep_times(ep+1);
end
save(fullfile(basepath,[basename,'.session.mat']),'session')
end



