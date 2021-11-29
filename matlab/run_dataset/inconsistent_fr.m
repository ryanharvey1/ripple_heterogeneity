% inconsistent_fr
% find inconsistent firing rate over time
% if firing rate is 0 for >50% of the recording session, tag as bad

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'Kenji'),:);
basepaths = unique(df.basepath);

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    
    start = cell_metrics.general.epochs{1}.startTime;
    stop = cell_metrics.general.epochs{end}.stopTime;
    bin_width = 60;
    bins = start:bin_width:stop;
    
    prop_zero = [];
    for cell_i = 1:length(cell_metrics.spikes.times)
        binned_st = histcounts(cell_metrics.spikes.times{cell_i},'BinEdges',bins) / bin_width;
        zerofrsum = sum(binned_st == 0);
        total_bins = length(binned_st);
        
        prop_zero(cell_i) = zerofrsum/total_bins;
    end
    if ~isfield(cell_metrics.tags,'Bad')
        cell_metrics.tags.Bad = [];
    end
    cell_metrics.tags.Bad = unique([cell_metrics.UID(prop_zero > .5),...
        cell_metrics.tags.Bad]);
    
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')

end