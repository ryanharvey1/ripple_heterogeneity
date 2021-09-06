% update_ab_data
df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepath = unique(df(contains(df.animal,'AB'),:).basepath);

for i = 1:length(basepath)
    disp(basepath{i})
    check_and_add(basepath{i})
end
for i = 1:length(basepath)
    disp(basepath{i})
    add_region(basepath{i})
end

function add_region(basepath)
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    try
        cell_metrics2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo_old.mat']));
    catch
        cell_metrics2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.legacy.mat']));
    end
    cell_metrics.brainRegion = cell_metrics2.cell_metrics.brainRegion;
    cell_metrics.CA1depth(1,:) = cell_metrics2.cell_metrics.CA1depth;
    
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end

function check_and_add(basepath)
basename = basenameFromBasepath(basepath);
% load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
load(fullfile(basepath,[basename,'.session.mat']))

pull_old_cell_metrics = false;
adjusted = false;

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
    adjusted = true;
end

if pull_old_cell_metrics
    cell_metrics_2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo1.mat']));
    spikes.filtWaveform = cell_metrics_2.cell_metrics.waveforms;
end

% check for nan waveforms
detected_nan = false;
for s = 1:length(spikes.UID)
    idx = isnan(spikes.filtWaveform{s});
    spikes.filtWaveform{s}(idx) = 0;
    if any(idx)
        detected_nan(s) = true;
    end
end
detected_nan = any(detected_nan);

if ~isfield(spikes,'peakVoltage') || detected_nan
    for s = 1:length(spikes.UID)
        spikes.peakVoltage(s) = max(spikes.filtWaveform{s}) -...
            min(spikes.filtWaveform{s});
    end
    adjusted = true;
end
if ~isfield(spikes,'timeWaveform') || detected_nan
    wf = [];
    for w = 1:length(spikes.filtWaveform)
        wf(:,w) = spikes.filtWaveform{w};
    end
    [~,midx] = min(wf);
    timeWaveform = ((1:size(wf,1)) - mode(midx)) / session.extracellular.sr * 1000;
    for s = 1:length(spikes.UID)
        spikes.timeWaveform{s} = timeWaveform;
    end
    adjusted = true;
end

if ~isfield(spikes,'channels_all')
    for s = 1:length(spikes.UID)
        spikes.channels_all{s} = session.extracellular.electrodeGroups.channels{spikes.shankID(s)};
    end
    adjusted = true;
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
    adjusted = true;
end

if ~isfield(spikes,'filtWaveform_all') || detected_nan
    for s = 1:length(spikes.UID)
        n_channels = session.extracellular.nChannels;
        spikes.filtWaveform_all{s} = zeros(n_channels,length(spikes.filtWaveform{s}));
        if isnan(spikes.maxWaveformCh1(s))
            spikes.filtWaveform_all{s}(1,:) = spikes.filtWaveform{s};
            spikes.filtWaveform_all{s} = imgaussfilt(spikes.filtWaveform_all{s},1);
            spikes.filtWaveform_all{s}(1,:) = spikes.filtWaveform{s};
        else
            spikes.filtWaveform_all{s}(spikes.maxWaveformCh1(s),:) = spikes.filtWaveform{s};
            spikes.filtWaveform_all{s} = imgaussfilt(spikes.filtWaveform_all{s},1);
            spikes.filtWaveform_all{s}(spikes.maxWaveformCh1(s),:) = spikes.filtWaveform{s};
        end
    end
    adjusted = true;
end

if ~isfield(spikes,'filtWaveform_all_std') || detected_nan
    for s = 1:length(spikes.UID)
        spikes.filtWaveform_all_std{s} = zeros(1,length(spikes.filtWaveform{s}));
    end
    adjusted = true;
end

% just in case... sort spikes. some sessions had some missorted spikes
for s = 1:length(spikes.UID)
    spikes.times{s} = sort(spikes.times{s});
end

if adjusted
    save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
end

if exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') && adjusted
    movefile(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),...
        fullfile(basepath,[basename,'.cell_metrics.cellinfo_old.mat']));
end

if adjusted || ~exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file')
    cell_metrics = ProcessCellMetrics('basepath',basepath,...
        'showGUI',false,...
        'spikes',spikes,...
        'getWaveformsFromDat',false,...
        'manualAdjustMonoSyn',false,...
        'session',session);
end
close all
end
