function run_ripple_pipe_girardeau(basepath,basename,spikes)

load(fullfile(basepath,[basename,'.session.mat']))

% check if emg exist
if ~exist(fullfile(basepath,[basename, '.EMGFromLFP.LFP.mat']),'file')
    bz_EMGFromLFP_km(basepath,...
        'samplingFrequency',10,'savemat',true,'noPrompts',true);
end

f_parts = strsplit(basepath,filesep);
load(fullfile(f_parts{1:end-1},'RecordingLocations.mat'))

idx = find(str2double(extractAfter(f_parts{end-1},'Rat')) == [shankInfo{2:end,end}])+1;
chan_to_check = [session.extracellular.electrodeGroups.channels{unique(shankInfo{idx,end-1}(:,2))}];

cell_idx = ismember(spikes.shankID,unique(shankInfo{idx,end-1}(:,2)));

event_threshold = true;
if ~any(cell_idx)
    disp('no hpc units')
    event_threshold = false;
end

% load lfp
lfp = getLFP(chan_to_check,'basepath',basepath,'basename',basename,'noPrompts',true);

EMGfilename = fullfile(basepath,[basename,'.EMGFromLFP.LFP.mat']);
if exist(EMGfilename)
    load(EMGfilename)
end

% estimate ripple band power normalized by wide band power
disp('finding ripple channel...')
try
    pBand = bandpower(double(lfp.data),lfp.samplingRate,[100,250]);
    pTot = bandpower(double(lfp.data),lfp.samplingRate,[1,(lfp.samplingRate/2)-1]);
catch
    for c = 1:size(lfp.data,2)
        pBand(c) = bandpower(double(lfp.data(:,c)),lfp.samplingRate,[100,250]);
        pTot(c) = bandpower(double(lfp.data(:,c)),lfp.samplingRate,[1,(lfp.samplingRate/2)-1]);
    end
end
[~,c_idx] = max(pBand./pTot);
c = chan_to_check(c_idx);

ripples = detect_ripple(lfp.data(:,c_idx),...
    lfp.timestamps,...
    restrict_spikes(spikes,cell_idx),...
    basepath,...
    EMGFromLFP,...
    event_threshold);

ripples.detectorinfo.detectionparms.ripple_channel = c;

save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end

function ripples = detect_ripple(lfp,timestamps,spikes,basepath,EMGFromLFP,event_threshold)

ripples = FindRipples(lfp,timestamps,...
    'durations',[50,300],...
    'thresholds',[1,3],...
    'passband',[100,250],...
    'minDuration',20,...
    'basepath',basepath,...
    'EMGFromLFP',EMGFromLFP);

if event_threshold
    ripples = eventSpikingTreshold(ripples,...
        'spikes',spikes,...
        'spikingThreshold',0.5,...
        'figOpt',false,...
        'basepath',basepath);
end
end

function temp_spikes = restrict_spikes(spikes,cell_idx)
temp_spikes = spikes;
n_cells = temp_spikes.numcells;

for f = fields(temp_spikes)'
    f = f{:};
    if [isnumeric(temp_spikes.(f)) ||...
            iscell(temp_spikes.(f))] &&...
            length(temp_spikes.(f)) == n_cells
        temp_spikes.(f) = temp_spikes.(f)(cell_idx);
    end
end
end