df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'OML'),:);

for i = 1:length(unique(df.basepath))
    basepath = df.basepath{i};
    disp(basepath)
    main(basepath)
end

function main(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.session.mat']))
% gui_session(session)
for ep = 1:length(session.epochs)
    if contains(session.epochs{ep}.name,'sleep','IgnoreCase',true)
        session.epochs{ep}.environment = 'sleep';
    elseif contains(session.epochs{ep}.name,'pre','IgnoreCase',true)
        session.epochs{ep}.environment = 'sleep';
    elseif contains(session.epochs{ep}.name,'post','IgnoreCase',true)
        session.epochs{ep}.environment = 'sleep';
    elseif contains(session.epochs{ep}.name,'baseline','IgnoreCase',true)
        session.epochs{ep}.environment = 'cheeseboard';
    elseif contains(session.epochs{ep}.name,'linear','IgnoreCase',true)
        session.epochs{ep}.environment = 'linear';
    elseif ~isfield(session.epochs{ep},'environment') 
        disp(session.epochs{ep}.name)
        session.epochs{ep}.environment = 'cheeseboard';
    elseif isempty(session.epochs{ep}.environment)
        disp(session.epochs{ep}.name)
        session.epochs{ep}.environment = 'cheeseboard';
    end
end
save(fullfile(basepath,[basename,'.session.mat']),'session')

end