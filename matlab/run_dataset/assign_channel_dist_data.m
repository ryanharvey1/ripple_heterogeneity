% assign_channel_dist_data

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\deep_sup_dist_estimation.csv');

idx = contains(df.polarity_reversal,'True');

plot(df.channelDistance(idx),df.channelDistance_pred(idx),'.k');

plot(df.channelDistance(~idx),df.channelDistance_pred(~idx),'.k');

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    
    if any(~contains(df(contains(df.basepath,basepath),:).polarity_reversal,'True'))
        load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']))
        temp_df = df(contains(df.basepath,basepath),:);
        
        idx = contains(temp_df.polarity_reversal,'True');
        figure;
        plot(temp_df.channelDistance(idx),temp_df.channelDistance_pred(idx),'.k');
        figure;

        plot(temp_df.channelDistance(~idx),temp_df.channelDistance_pred(~idx),'.k');
    end

end