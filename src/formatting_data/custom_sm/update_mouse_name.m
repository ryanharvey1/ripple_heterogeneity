% update_mouse_name
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\mouse_sessions.csv');
for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    animName = animalFromBasepath(basepath);
    session.animal.name = animName;
    
    % update species and strain while we are here
    session.animal.species = 'Mouse';
    session.animal.strain = 'C57B1/6';
    
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end    
