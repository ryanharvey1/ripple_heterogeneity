% run_rank_order

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);


WaitMessage = parfor_wait(length(basepaths));
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    main(basepath,basename)
    WaitMessage.Send;
end
WaitMessage.Destroy;

function main(basepath,basename)

if exist(fullfile(basepath,[basename,'.rankStats.mat']),'file')
    load(fullfile(basepath,[basename,'.rankStats.mat']))
    if sum(isnan(rankStats.rankClusters)) == 0
        return
    end
end

load(fullfile(basepath,[basename,'.ripples.events.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

ripSpk = bz_getRipSpikes('basepath',basepath,...
                        'basename',basename,...
                        'spikes',spikes,...
                        'events',ripples.timestamps,...
                        'saveMat',false); 
                    
[rankStats] = RankOrder('basepath',basepath,...
    'spkEventTimes',ripSpk);
end