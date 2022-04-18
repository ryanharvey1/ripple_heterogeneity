% fix_UID_bug

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);
basepaths = flipud(basepaths);

incorrect_UID = [];
for i = 1:length(basepaths)
    basepath = basepaths{i};
    disp(basepath)
    basename = basenameFromBasepath(basepath);
    
    % check to see if UID increases by 1 step
    % cluID should be the original
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']));
    
    if sum(diff(spikes.UID) == 1) ~= spikes.numcells-1
        
        incorrect_UID = [incorrect_UID;{basepath}];
        
        spikes.cluID = spikes.UID;
        spikes.UID = 1:spikes.numcells;
        spikes = get_spindices(spikes);
        
        save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes');

        ProcessCellMetrics('basepath',basepath,...
            'spikes',spikes,...
            'manualAdjustMonoSyn',false,...
            'showGUI',false)
        
        close all
    end
end

