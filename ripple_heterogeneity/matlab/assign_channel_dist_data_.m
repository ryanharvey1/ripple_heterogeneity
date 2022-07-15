% assign_channel_dist_data_

% df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\deep_sup_dist_estimation.csv');
% df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\deep_sup_dist_estimation_v2_6_10_22.csv');
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\deep_sup_dist_estimation_v2_7_15_22.csv');

basepaths = unique(df.basepath);
for i = 1:length(unique(df.basepath))
    basepath = basepaths{i};
    % detect if some shanks don't have a polarity reversal
    if any(~contains(df(contains(df.basepath,basepath),:).polarity_reversal,'True'))
        % reassign pyr distance
        disp(basepath)
        add_new_distance(df,basepath)
    end
end

function add_new_distance(df,basepath)
basename = basenameFromBasepath(basepath);

% need chan coords to run
load(fullfile(basepath,[basename,'.session.mat']),'session')
if ~isfield(session.extracellular,'chanCoords')
    return
end
if isempty(session.extracellular.chanCoords.x)
   return 
end

load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'deepSuperficialfromRipple')

% check if already corrected
if isfield(deepSuperficialfromRipple.processinginfo.params,'adjusted_by_predictive_model')
    if deepSuperficialfromRipple.processinginfo.params.adjusted_by_predictive_model
        return 
    end
end

load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')

cell_metrics.general.SWR = deepSuperficialfromRipple;
deepSuperficial_ChClass = deepSuperficialfromRipple.channelClass;
cell_metrics.general.deepSuperficial_file =...
    fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']);

% save backup
if ~exist(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo_old_61022.mat']),'file')
    save(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo_old_61022.mat']),...
        'deepSuperficialfromRipple')
end

temp_df = df(contains(df.basepath,basepath),:);
new_df = table();
for i = unique(temp_df.shank)'
    cur_df = temp_df(temp_df.shank == i,:);
    [a,b] = ismember(deepSuperficialfromRipple.ripple_channels{i},cur_df.channel);
    new_df = [new_df;cur_df(b,:)];
end

% vert_space = deepSuperficialfromRipple.processinginfo.params.verticalSpacing;
for i = unique(new_df.shank)'
    
    channelClass = unique(new_df(new_df.shank == i,:).channelClass);
    if length(channelClass) == 1
        if contains(channelClass,"Superficial")
            new_dist = new_df(new_df.shank == i,:).channelDistance_pred;   
            % offset by first estimated channel distance
            y = session.extracellular.chanCoords.y(new_df(new_df.shank == i,:).channel);
            y = abs(y - max(y));
            new_df(new_df.shank == i,:).channelDistance = y + new_dist(1);
            
        elseif contains(channelClass,"Deep")
            new_dist = new_df(new_df.shank == i,:).channelDistance_pred;
            % offset by first estimated channel distance
            y = session.extracellular.chanCoords.y(new_df(new_df.shank == i,:).channel);
            y = abs(y - max(y));
            new_dist = y - max(y) + new_dist(1);
            new_df(new_df.shank == i,:).channelDistance = new_dist;
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

deepSuperficialfromRipple.processinginfo.params.adjusted_by_predictive_model = true;

save(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),...
    'deepSuperficialfromRipple')
save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end
