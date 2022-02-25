% add_ripples

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\mouse_sessions.csv');
for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    ripple_file = fullfile(basepath,[basename,'.ripples.events.mat']);
    if exist(ripple_file,'file')
        continue
    end
    cd(basepath)
    % get ripple params from previous file and detect again
    ripple_file = fullfile(basepath,[basename,'.SWR.events.mat']);
    load(ripple_file)
    
    SWR.params
end