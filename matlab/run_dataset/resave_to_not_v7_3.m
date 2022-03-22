% resave_to_not_v7_3

df = readtable("Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv")
for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath, [basename,'.deepSuperficialfromRipple.channelinfo.mat']));
    save(fullfile(basepath, [basename,'.deepSuperficialfromRipple.channelinfo.mat']),'deepSuperficialfromRipple');
end
