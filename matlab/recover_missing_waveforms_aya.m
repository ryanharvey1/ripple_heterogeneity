% recover_missing_waveforms_aya


df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepath = unique(df(contains(df.animal,'AYA'),:).basepath);

for i = 1:length(basepath)
    if exist(fullfile(basepath{i},[basenameFromBasepath(basepath{i}),'.dat']),'file')
        dat_exist(i) = true;
    else
        dat_exist(i) = false;
    end
end

basepath = basepath(dat_exist);

for i = 1:length(basepath)
    disp(basepath{i})
    try
        load(fullfile(basepath{i},[basenameFromBasepath(basepath{i}),'.spikes.cellinfo.mat']))
        load(fullfile(basepath{i},[basenameFromBasepath(basepath{i}),'.cell_metrics.cellinfo.mat']))
        load(fullfile(basepath{i},[basenameFromBasepath(basepath{i}) '.session.mat']));
    catch
        continue
    end
    % check for nan waveforms
    extract_from_dat = false;
    for w = 1:length(cell_metrics.waveforms.filt)
        if any(isnan(cell_metrics.waveforms.filt{1, w}))
            extract_from_dat = true;
            break
        end
    end
    if extract_from_dat
        disp(['extracting waveforms from: ', basepath{i}])
        
        spikes = getWaveformsFromDat(spikes,session);
        save(fullfile(basepath{i},[basenameFromBasepath(basepath{i}),'.spikes.cellinfo.mat']),'spikes')
        
        % rename cell metrics so we can load from scratch
        movefile(fullfile(basepath{i},[basenameFromBasepath(basepath{i}),'.cell_metrics.cellinfo.mat']),...
            fullfile(basepath{i},[basenameFromBasepath(basepath{i}),'.cell_metrics.cellinfo_old.mat']));
        
        cell_metrics = ProcessCellMetrics('basepath',basepath{i},...
            'showGUI',false,...
            'spikes',spikes,...
            'getWaveformsFromDat',true,...
            'manualAdjustMonoSyn',false,...
            'session',session);
    end
end
