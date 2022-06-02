% custom_for_fujisawaS
% folders = dir('Z:\Data\FujisawaS\EE')

basepaths = {'Z:\Data\FujisawaS\EE\EE0622fm',...
    'Z:\Data\FujisawaS\EE\EE0627fm',...
    'Z:\Data\FujisawaS\EE\EE0705fm',...
    'Z:\Data\FujisawaS\EE\EE0706fm',...
    'Z:\Data\FujisawaS\EE\EE0708fm',...
    }

%% concat lfp files
for i = 1:length(basepaths)
    dat_files = dir(fullfile(basepaths{i},'Subsessions','**','*.dat'));
    eeg_files = dir(fullfile(basepaths{i},'Subsessions','**','*.eeg'));
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    newdatpath = fullfile(basepath,[basename,'.lfp']);
    if ~exist(newdatpath,'file')
        datpathsplus = [];
        for didx = 1:length(eeg_files)
            if didx == length(eeg_files)
                datpathsplus{didx} = [fullfile(eeg_files(didx).folder,eeg_files(didx).name)];
            else
                datpathsplus{didx} = [fullfile(eeg_files(didx).folder,eeg_files(didx).name) ' +'];
            end
        end
        
        cs = strjoin(datpathsplus);
        catstring = ['! copy /b ', cs, ' ',newdatpath];
        
        disp('Concatenating Amplifier Dats... be patient')
        eval(catstring) %execute concatention
        
        t = dir(newdatpath);
        if t.bytes ~= sum(vertcat(eeg_files.bytes))
            error('New .dat size not right.  Exiting')
        end
    end
end
%% make merge points
for i = 1:length(basepaths)
    dat_files = dir(fullfile(basepaths{i},'Subsessions','**','*.dat'));
    eeg_files = dir(fullfile(basepaths{i},'Subsessions','**','*.eeg'));

    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    sessionInfo = LoadXml(fullfile(basepath,[basename, '.xml']));
    nSamp = [];
    for didx = 1:length(eeg_files)
        % use dir to get size of dat file in bytes
        % determine number of bytes per sample
        dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); 
        % Number of samples per channel
        nSamp(didx) = eeg_files(didx).bytes/(sessionInfo.nChannels*dataTypeNBytes);  
    end
    cumsum_nSamp = cumsum([nSamp]);
    starts = [0,cumsum_nSamp(1:end-1)];
    transitiontimes_samp = [starts',cumsum_nSamp'];
    transitiontimes_sec = transitiontimes_samp./sessionInfo.lfpSampleRate;
    
    firstlasttimepoints = [zeros(length(nSamp),1),nSamp'];
    
    recordingnames = [];
    for didx = 1:length(eeg_files)
        recordingnames{1,didx} = basenameFromBasepath(eeg_files(didx).folder);
    end
    
    MergePoints.timestamps = transitiontimes_sec;
    MergePoints.timestamps_samples = transitiontimes_samp;
    MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
    MergePoints.foldernames = recordingnames;
    MergePoints.detectorinfo.detectorname = 'custom_for_fujisawaS';
    MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
    
    save(fullfile(basepath,[basename,'.MergePoints.events.mat']),'MergePoints');
end


%% make basename.session
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    session = sessionTemplate(basepath,'showGUI',true);
    save(fullfile(basepath,[basename, '.session.mat']),'session');
end

%% sleep states
for i = 1:length(basepaths)
    basepath = basepaths{i};

    SleepScoreMaster(basepath,'noPrompts',true); % takes lfp in base 0
    thetaEpochs(basepath);
end