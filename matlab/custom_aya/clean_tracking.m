df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df = df(contains(df.basepath,'AYAold'),:);
% df = df(contains(df.basepath,'GrosmarkAD'),:);
df = df(~contains(df.basepath,'GirardeauG') &...
            ~contains(df.basepath,'Kenji'),:);


for i = 1:length(unique(df.basepath))
    basepath = df.basepath{i};
    disp(basepath)
    main(basepath)
end

function main(basepath)
basename = basenameFromBasepath(basepath);

if exist(fullfile(basepath,[basename,'.restrictxy.mat']),'file')
    disp(['Finished'])
    
    load(fullfile(basepath,[basename,'.restrictxy.mat']))
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))

    behavior.position.x(~good_idx) = NaN;
    behavior.position.y(~good_idx) = NaN;
    
    save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
    return
end

load(fullfile(basepath,[basename,'.session.mat']))
load(fullfile(basepath,[basename,'.animal.behavior.mat']))

if ~isempty(behavior.position.x)
    start = [];
    stop = [];
    
    for ep = 1:length(session.epochs)
        if ~contains(session.epochs{ep}.environment,'sleep')
            start = [start,session.epochs{ep}.startTime];
            stop = [stop,session.epochs{ep}.stopTime];
        end
    end
    
    good_idx = manual_trackerjumps(behavior.time,...
        behavior.position.x,...
        behavior.position.y,...
        start,...
        stop,...
        basepath);
    
    behavior.position.x(~good_idx) = NaN;
    behavior.position.y(~good_idx) = NaN;
    
    save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end
% for ep = 1:length(session.epochs)
%     if ~contains(session.epochs{ep}.environment,'sleep')
%         disp(session.epochs{ep}.environment)
%         start = session.epochs{ep}.startTime;
%         stop = session.epochs{ep}.stopTime;
%         idx = behavior.time >= start & behavior.time <= stop;
%         figure;
%         plot(behavior.position.x(idx),behavior.position.y(idx))
%     end
% end
end



