
% correct_session_basepaths
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    load(fullfile(basepath,[basename,'.session.mat']))
    session.general.basePath = basepath;
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end
