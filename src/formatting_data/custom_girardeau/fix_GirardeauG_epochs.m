% fix_GirardeauG_epochs
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'GirardeauG'),:);
basepaths = unique(df.basepath);
parfor i = 1:length(basepaths)
    basepath = basepaths{i};
    disp(basepath)
    main(basepath)
end

function main(basepath)

basename = basenameFromBasepath(basepath);

if ~exist(fullfile(basepath,[basename,'.cat.evt']),'file')
    return
end
load(fullfile(basepath,[basename,'.session.mat']));

events = LoadEvents(fullfile(basepath,[basename,'.cat.evt']));

% iter through .cat.evt and use those epochs to fill .session
session.epochs = [];
ep = 1;
for e = 1:2:(length(events.description))
    times = events.time(e:e+1);
    session.epochs{ep}.startTime = times(1);
    session.epochs{ep}.stopTime = times(2);
    name = strsplit(events.description{e},'-');
    session.epochs{ep}.name = name{end};
    disp(session.epochs{ep}.name)
    ep = ep+1;
end

load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']));

cell_metrics = ProcessCellMetrics('basepath',basepath,...
    'showGUI',false,...
    'spikes',spikes,...
    'getWaveformsFromDat',false,...
    'manualAdjustMonoSyn',false,...
    'session',session);

save(fullfile(basepath,[basename,'.session.mat']),'session');
save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics');

close all
end