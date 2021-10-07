% update_cell_metrics_deep_superficial
%
% checks to see if cell_metrics.deepSuperficial has been updated by deepSuperficialfromRipple.channelinfo
% if it has not been updated...this will update it.

df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepaths = unique(df.basepath);

overwrite = true;

for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    disp(basepath)
    
    if all(contains(cell_metrics.deepSuperficial,"Unknown")) || overwrite
        
        deepSuperficial_file = fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']);
        
        load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
        load(deepSuperficial_file,'deepSuperficialfromRipple')
        
        cell_metrics.general.SWR = deepSuperficialfromRipple;
        deepSuperficial_ChDistance = deepSuperficialfromRipple.channelDistance;
        deepSuperficial_ChClass = deepSuperficialfromRipple.channelClass;
        cell_metrics.general.deepSuperficial_file = deepSuperficial_file;
        
        for j = 1:cell_metrics.general.cellCount
            if isnan(spikes.maxWaveformCh1(j))
                cell_metrics.deepSuperficial(j) = {'Unknown'};
                cell_metrics.deepSuperficialDistance(j) = NaN;
            else
                cell_metrics.deepSuperficial(j) =...
                    deepSuperficial_ChClass(spikes.maxWaveformCh1(j)); % cell_deep_superficial OK
                cell_metrics.deepSuperficialDistance(j) =...
                    deepSuperficial_ChDistance(spikes.maxWaveformCh1(j)); % cell_deep_superficial_distance
            end
        end
        save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    end
end

