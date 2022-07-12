
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
transitiontimes_sec = transitiontimes_samp./sessionInfo.lfpSampleRate;

firstlasttimepoints = [zeros(length(nSamp),1),nSamp'];

recordingnames = [];
for didx = 1:length(dat_files)
    recordingnames{1,didx} = basenameFromBasepath(dat_files(didx).folder);
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

%% get cell metrics

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,...
        'manualAdjustMonoSyn',false,'showGUI',false,'getWaveformsFromDat',false);
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