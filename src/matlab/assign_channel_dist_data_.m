% assign_channel_dist_data_

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\deep_sup_dist_estimation.csv');

basepaths = unique(df.basepath);
parfor i = 1:length(unique(df.basepath))
    basepath = basepaths{i};
    % detect if some shanks don't have a polarity reversal
    if any(~contains(df(contains(df.basepath,basepath),:).polarity_reversal,'True'))
        % reassign pyr distance
        add_new_distance(df,basepath)
    end
end

function add_new_distance(df,basepath)
basename = basenameFromBasepath(basepath);

load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))

cell_metrics.general.SWR = deepSuperficialfromRipple;
deepSuperficial_ChClass = deepSuperficialfromRipple.channelClass;
cell_metrics.general.deepSuperficial_file =...
    fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']);

% save backup
if ~exist(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo_old.mat']),'file')
    save(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo_old.mat']),...
        'deepSuperficialfromRipple')
end

temp_df = df(contains(df.basepath,basepath),:);
new_df = table();
for i = unique(temp_df.shank)'
    cur_df = temp_df(temp_df.shank == i,:);
    [a,b] = ismember(deepSuperficialfromRipple.ripple_channels{i},cur_df.channel);
    new_df = [new_df;cur_df(b,:)];
end

vert_space = deepSuperficialfromRipple.processinginfo.params.verticalSpacing;
for i = unique(new_df.shank)'
    
    channelClass = unique(new_df(new_df.shank == i,:).channelClass);
    if length(channelClass) == 1
        if contains(channelClass,"Superficial")
            new_dist = new_df(new_df.shank == i,:).channelDistance_pred;
            new_dist = linspace(0,length(new_dist)*vert_space,length(new_dist)) + new_dist(1);
            new_df(new_df.shank == i,:).channelDistance = new_dist';
        elseif contains(channelClass,"Deep")
            new_dist = new_df(new_df.shank == i,:).channelDistance_pred;
            new_dist = linspace(-length(new_dist)*vert_space,0,length(new_dist)) + new_dist(1);
            new_df(new_df.shank == i,:).channelDistance = new_dist';
        end
    end
end

[~,Lib] = ismember(new_df.channel,deepSuperficialfromRipple.channel);
deepSuperficialfromRipple.channelDistance(Lib) = new_df.channelDistance';

for j = 1:cell_metrics.general.cellCount
    if isnan(spikes.maxWaveformCh1(j))
        cell_metrics.deepSuperficial(j) = {'Unknown'};
        cell_metrics.deepSuperficialDistance(j) = NaN;
    else
        cell_metrics.deepSuperficial(j) =...
            deepSuperficial_ChClass(spikes.maxWaveformCh1(j));
        cell_metrics.deepSuperficialDistance(j) =...
            deepSuperficialfromRipple.channelDistance(spikes.maxWaveformCh1(j));
    end
end

save(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),...
    'deepSuperficialfromRipple')
save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end
