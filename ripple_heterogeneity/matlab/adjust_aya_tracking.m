% adjust_aya_tracking

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'AYAold'),:);

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
    if ~isempty(behavior.position.x)
        figure;
        plot(behavior.position.x,behavior.position.y)
        title(basepath)
        if input('re run? ') == 1
            behavior = general_behavior_file('basepath',basepath,'save_mat',false);
            figure;
            plot(behavior.position.x,behavior.position.y,'g')
            title(basepath)
            if input('looks better? ')
                save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
            end
        end
        close all
    end
end


