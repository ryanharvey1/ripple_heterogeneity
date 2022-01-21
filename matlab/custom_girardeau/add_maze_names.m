df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'GirardeauG'),:);

for i = 1:length(unique(df.basepath))
    basepath = df.basepath{i};
    disp(basepath)
    main(basepath)
end

function main(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.session.mat']))

for ep = 1:length(session.epochs)
    if contains(session.epochs{ep}.name,'sleep')
        session.epochs{ep}.environment = 'sleep';
    elseif contains(session.epochs{ep}.name,'run')
        session.epochs{ep}.environment = 'linear';
    elseif contains(session.epochs{ep}.name,'water')
        session.epochs{ep}.environment = 'water';
    elseif contains(session.epochs{ep}.name,'smell')
        session.epochs{ep}.environment = 'smell';
    elseif contains(session.epochs{ep}.name,'video')
        session.epochs{ep}.environment = 'unknown';
    else
        keyboard
    end
end
save(fullfile(basepath,[basename,'.session.mat']),'session')

end