% inconsistent_fr
% find inconsistent firing rate over time
% if firing rate is 0 for >40% of the recording session, tag as bad

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'Kenji'),:);
basepaths = unique(df.basepath);

prop_zero_thres = .4;
check_epoch_fr = true;

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    
    start = cell_metrics.general.epochs{1}.startTime;
    stop = cell_metrics.general.epochs{end}.stopTime;
    bin_width = 60;
    bins = start:bin_width:stop;
    
    % check for total percent that is zero fr
    prop_zero = [];
    for cell_i = 1:length(cell_metrics.spikes.times)
        binned_st = histcounts(cell_metrics.spikes.times{cell_i},'BinEdges',bins) / bin_width;
        zerofrsum = sum(binned_st == 0);
        total_bins = length(binned_st);
        
        prop_zero(cell_i) = zerofrsum/total_bins;
    end
    if ~isfield(cell_metrics,'tags')
        cell_metrics.tags.Bad = [];
    end
    if ~isfield(cell_metrics.tags,'Bad')
        cell_metrics.tags.Bad = [];
    end
    cell_metrics.tags.Bad = unique([cell_metrics.UID(prop_zero > prop_zero_thres),...
        cell_metrics.tags.Bad]);
    
   
    % make sure fr doesn't drop down to zero for an epoch
    prop_zero_cell = [];
    bin_width = 15;
    if check_epoch_fr
        for cell_i = 1:length(cell_metrics.spikes.times)
            prop_zero = [];
            for ep = 1:length(cell_metrics.general.epochs)
                st = cell_metrics.spikes.times{cell_i};
                idx = st >= cell_metrics.general.epochs{ep}.startTime &...
                    st <= cell_metrics.general.epochs{ep}.stopTime;
                start = cell_metrics.general.epochs{ep}.startTime;
                stop = cell_metrics.general.epochs{ep}.stopTime;
                bins = start:bin_width:stop;
                binned_st = histcounts(st(idx),'BinEdges',bins) / bin_width;
                zerofrsum = sum(binned_st == 0);
                total_bins = length(binned_st);
                prop_zero(ep) = zerofrsum/total_bins;
            end
            prop_zero_cell(cell_i) = any(prop_zero==1);
        end
    end
    cell_metrics.tags.Bad = unique([cell_metrics.UID(logical(prop_zero_cell)),...
        cell_metrics.tags.Bad]);
    
    
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    
    total_units(i) = length(cell_metrics.UID);
    bad_unit_count(i) = length(cell_metrics.tags.Bad);
end