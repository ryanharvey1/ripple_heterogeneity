
basepaths = {'Z:\Data\ORproject\OR15\hc280118','Z:\Data\ORproject\OR15\hc300118'}


%% make merge points
% for i = 1:length(basepaths)
basepath = pwd;

dat_files = dir(fullfile(basepath,'**','*amplifier.dat'));
% eeg_files = dir(fullfile(basepath,'**','*.eeg'));


basename = basenameFromBasepath(basepath);

sessionInfo = LoadXml(fullfile(basepath,[basename, '.xml']));
nSamp = [];
for didx = 1:length(dat_files)
    % use dir to get size of dat file in bytes
    % determine number of bytes per sample
    dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8'));
    % Number of samples per channel
    nSamp(didx) = dat_files(didx).bytes/(sessionInfo.nChannels*dataTypeNBytes);
end
cumsum_nSamp = cumsum([nSamp]);
starts = [0,cumsum_nSamp(1:end-1)];
transitiontimes_samp = [starts',cumsum_nSamp'];
transitiontimes_sec = transitiontimes_samp./sessionInfo.SampleRate;

firstlasttimepoints = [zeros(length(nSamp),1),nSamp'];

recordingnames = [];
for didx = 1:length(dat_files)
    [a,b] = fileparts(dat_files(didx).folder);
    recordingnames{1,didx} = b;
end

MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.detectorinfo.detectorname = 'custom_for_fujisawaS';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

save(fullfile(basepath,[basename,'.MergePoints.events.mat']),'MergePoints');
% end

%%
session = sessionTemplate(pwd,'showGUI',true);

coords = readtable("D:\github\ripple_heterogeneity\data\electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177 (1).csv");

chanCoords.x(1:64) = coords.x_um_;
chanCoords.y(1:64) = coords.y_um_;

%% sleep states
for i = 1:length(basepaths)
    basepath = basepaths{i};

    SleepScoreMaster(basepath,'noPrompts',true); % takes lfp in base 0
    thetaEpochs(basepath);
end

%% get spikes
for i = 1:length(basepaths)
    basepath = basepaths{i};
    spikes = loadSpikes('basepath',basepath)
%     spikes = loadSpikes('basepath',basepath,'getWaveformsFromSource',true,'getWaveformsFromDat',false);
end

%% ripples
basepath = pwd
basename = basenameFromBasepath(basepath);

load(fullfile(basepath,[basename,'.session.mat']))
ripples = DetectSWR([session.channelTags.Ripple.channels, session.channelTags.SharpWave.channels],...
    'basepath',basepath,...
    'saveMat',true,'thresSDswD', [0.25, 1],'thresSDrip', [0.25, 1],...
    'forceDetect',true,'check',true);


%% opto pulse
old_pulses = load('hc300118.stimPulses.mat');
% old_pulses = load('hc280118.stimPulses.mat');
% basepaths = {'Z:\Data\ORproject\OR15\hc280118','Z:\Data\ORproject\OR15\hc300118'}

% added closed loop and delayed pulses, pulses lasted for 100ms
pulses.timestamps(:,1) = [old_pulses.pulsesCL;old_pulses.pulsesD];
pulses.timestamps(:,2) = pulses.timestamps(:,1) + 0.1;

% add event label (closed loop is 0 and delayed is 1)
pulses.eventGroupID = [zeros(length(old_pulses.pulsesCL),1) ;...
    zeros(length(old_pulses.pulsesD),1) + 1];

% make sure ts are sorted
[~,idx] = sort(pulses.timestamps(:,1));
pulses.timestamps = pulses.timestamps(idx,:);


% pulses.timestamps(:,1) = old_pulses.pulses;
% pulses.timestamps(:,2) = old_pulses.pulses+.1;

pulses.amplitude = nan(length(pulses.timestamps),1);

% pulses.eventGroupID(ismember(pulses.timestamps(:,1),old_pulses.pulsesD)) = 0;
% pulses.eventGroupID(ismember(pulses.timestamps(:,1),old_pulses.pulsesCL)) = 1;

pulses.duration = pulses.timestamps(:,2) - pulses.timestamps(:,1);

optoStim.timestamps = pulses.timestamps;
optoStim.peaks = median(pulses.timestamps,2);
optoStim.amplitude = pulses.amplitude;
optoStim.amplitudeUnits = 'au';
optoStim.eventID = pulses.eventGroupID;
optoStim.eventIDlabels = {'closed_loop','delayed'};
optoStim.eventIDbinary = true;
optoStim.center = median(pulses.timestamps,2);
optoStim.duration = pulses.duration;
optoStim.detectorinfo = 'getAnalogPulses_pulses_to_manipulation';

saveStruct(optoStim,'manipulation','session',session);

%% get cell metrics

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
    
    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,...
        'manualAdjustMonoSyn',false,'showGUI',true);
    
%     cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,...
%         'manualAdjustMonoSyn',false,'showGUI',false,'getWaveformsFromDat',false);
end

% basepath = pwd
% basename = basenameFromBasepath(basepath);
% 
% load(fullfile(basepath,[basename,'.session.mat']))
% 
% ripples = DetectSWR([session.channelTags.Ripple.channels, session.channelTags.SharpWave.channels],...
%     'basepath',basepath,...
%     'saveMat',true,'thresSDswD', [0.25, 1],'thresSDrip', [0.25, 1],...
%     'forceDetect',true,'check',true);