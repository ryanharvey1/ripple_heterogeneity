function run_ripple_pipe_kenji(basepath,basename,spikes)

% load electrode position info
load('A:\Data\Kenji\ElePosition.mat')
shank_region = ElePosition(contains(ElePosition(:,2),basename),6:end);
for i = 1:length(shank_region);shank_region{i}=lower(shank_region{i});end

load(fullfile(basepath,[basename,'.session.mat']))

% check if emg exist
if ~exist(fullfile(basepath,[basename, '.EMGFromLFP.LFP.mat']),'file')
    bz_EMGFromLFP(basepath,'basename',basename,...
        'samplingFrequency',10,'savemat',true,'noPrompts',true);
end

% if data has ca1 then we can look at ripples that co-occur with sharp waves
if any(strcmp(shank_region,'ca1') | strcmp(shank_region,'ca1c'))
    
    % find out what channel has ca1 and run that channel
    idx = find(strcmp(shank_region,'ca1') | strcmp(shank_region,'ca1c'));
    
elseif any(strcmp(shank_region,'ca3') |...
        strcmp(shank_region,'ca2') |...
        strcmp(shank_region,'ca') |...
        strcmp(shank_region,'dg') |...
        strcmp(shank_region,'dgca3'))
    
    % find out what channel has HPC and run that channel
    idx = find(strcmp(shank_region,'ca3') |...
        strcmp(shank_region,'ca2') |...
        strcmp(shank_region,'ca') |...
        strcmp(shank_region,'dg') |...
        strcmp(shank_region,'dgca3'));
else
    return 
end

chan_to_check = [session.extracellular.electrodeGroups.channels{idx}];

% load lfp
disp('loading lfp...')
lfp = bz_GetLFP(chan_to_check-1,'basepath',basepath,'basename',basename,'noPrompts',true);

cell_idx = ismember(spikes.shankID,idx);

EMGfilename = fullfile(basepath,[basename,'.EMGFromLFP.LFP.mat']);
if exist(EMGfilename)
    load(EMGfilename)
end

% estimate ripple band power normalized by wide band power
disp('finding ripple channel...')

pBand = bandpower(double(lfp.data),lfp.samplingRate,[100,250]);
pTot = bandpower(double(lfp.data),lfp.samplingRate,[1,(lfp.samplingRate/2)-1]);
[~,c_idx] = max(pBand./pTot);
c = chan_to_check(c_idx);

ripples = detect_ripple(lfp.data(:,c_idx),...
    lfp.timestamps,...
    restrict_spikes(spikes,cell_idx),...
    basepath,...
    EMGFromLFP);

ripples.detectorinfo.detectionparms.ripple_channel = c;

save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end

function ripples = detect_ripple(lfp,timestamps,spikes,basepath,EMGFromLFP)

ripples = bz_FindRipples(lfp,timestamps,...
    'durations',[50,300],...
    'thresholds',[1,3],...
    'passband',[100,250],...
    'minDuration',20,...
    'basepath',basepath,...
    'EMGFromLFP',EMGFromLFP);

ripples = eventSpikingTreshold(ripples,...
    'spikes',spikes,...
    'spikingThreshold',0.5,...
    'figOpt',false,...
    'basepath',basepath);
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