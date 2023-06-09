%% make merge points
for i = 1:length(basepaths)
    dat_files = dir(fullfile(basepaths{i},'Subsessions','**','*.dat'));
    % eeg_files = dir(fullfile(basepaths{i},'Subsessions','**','*.eeg'));

    basepath = basepaths{i};
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
    cumsum_nSamp = cumsum(nSamp);
    starts = [0,cumsum_nSamp(1:end-1)];
    transitiontimes_samp = [starts',cumsum_nSamp'];
    transitiontimes_sec = transitiontimes_samp./sessionInfo.SampleRate;
   
    firstlasttimepoints = [zeros(length(nSamp),1),nSamp'];
    
    recordingnames = [];
    for didx = 1:length(dat_files)
        recordingnames{1,didx} = basenameFromBasepath(dat_files(didx).folder);
    end
    
    MergePoints.timestamps = transitiontimes_sec;
    MergePoints.timestamps_samples = transitiontimes_samp;
    MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
    MergePoints.foldernames = recordingnames;
    MergePoints.detectorinfo.detectorname = 'custom';
    MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
    
    save(fullfile(basepath,[basename,'.MergePoints.events.mat']),'MergePoints');
end