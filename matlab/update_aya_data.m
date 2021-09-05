% update_aya_data
df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepath = unique(df(contains(df.animal,'AYA'),:).basepath);

for i = 1:length(basepath)
    disp(basepath{i})
    check_and_add(basepath{i})
end

function check_and_add(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
load(fullfile(basepath,[basename,'.session.mat']))

pull_old_cell_metrics = false;

if ~isfield(spikes,'filtWaveform')
    load(fullfile(basepath,['cell_waveforms.mat']))
    for s = 1:length(spikes.UID)
        if s > size(waveforms,1)
            warning('size mix-up...pull from older cell_metrics')
            pull_old_cell_metrics = true;
            break
        end
        spikes.filtWaveform{s} = waveforms(s,:);
    end
end
if pull_old_cell_metrics
    cell_metrics_2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo1.mat']));
    spikes.filtWaveform = cell_metrics_2.cell_metrics.waveforms;
end
if ~isfield(spikes,'peakVoltage')
    for s = 1:length(spikes.UID)
        spikes.peakVoltage(s) = max(spikes.filtWaveform{s}) -...
            min(spikes.filtWaveform{s});
    end
end
if ~isfield(spikes,'timeWaveform')
    wf = [];
    for w = 1:length(spikes.filtWaveform)
        wf(:,w) = spikes.filtWaveform{w};
    end
    [~,midx] = min(wf);
    timeWaveform = ((1:size(waveforms,2)) - mode(midx)) / session.extracellular.sr * 1000;
    for s = 1:length(spikes.UID)
        spikes.timeWaveform{s} = timeWaveform;
    end
end
if ~isfield(spikes,'channels_all')
    for shank = 1:length(session.extracellular.electrodeGroups.channels)
        for s = 1:length(spikes.UID)
            spikes.channels_all{s} = session.extracellular.electrodeGroups.channels{shank};
        end
    end
end

if ~isfield(spikes,'maxWaveformCh') || ~isfield(spikes,'maxWaveformCh1')
    if pull_old_cell_metrics
        spikes.maxWaveformCh1 = cell_metrics_2.cell_metrics.maxWaveformCh1;
        spikes.maxWaveformCh = spikes.maxWaveformCh1-1;
    else
        load(fullfile(basepath,['cell_waveforms.mat']))
        for s = 1:length(spikes.UID)
            spikes.maxWaveformCh1(s) =...
                session.extracellular.electrodeGroups.channels{location(s,1)}(location(s,2));
            spikes.maxWaveformCh(s) = spikes.maxWaveformCh1(s)-1;
        end
    end
end
save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')

cell_metrics = ProcessCellMetrics('basepath',basepath,...
    'showGUI',false,...
    'spikes',spikes,...
    'getWaveformsFromDat',false,...
    'manualAdjustMonoSyn',false,...
    'session',session);

close all
end
