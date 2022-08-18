% custom_for_fujisawaS

basepaths = {'Z:\Data\FujisawaS\FF\FF1114',...
    'Z:\Data\FujisawaS\FF\FF1116',...
    'Z:\Data\FujisawaS\FF\FF1119',...
    'Z:\Data\FujisawaS\EE\EE0710',...
    'Z:\Data\FujisawaS\EE\EE0711',...
    'Z:\Data\FujisawaS\GG\GG0401',...
    'Z:\Data\FujisawaS\GG\GG0406'}

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

%% get spikes
basepaths = {'Z:\Data\FujisawaS\FF\FF1119',...
    'Z:\Data\FujisawaS\EE\EE0710',...
    'Z:\Data\FujisawaS\EE\EE0711',...
    'Z:\Data\FujisawaS\GG\GG0401',...
    'Z:\Data\FujisawaS\GG\GG0406'}
for i = 1:length(basepaths)
    basepath = basepaths{i};

    spikes = loadSpikes('basepath',basepath,...
        'getWaveformsFromSource',true,...
        'getWaveformsFromDat',true,...
        'format','klustakwik','forceReload',true);
end


%% get cell metrics

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,...
        'manualAdjustMonoSyn',false,'showGUI',true,...
        'excludeManipulationIntervals',false);
end

%% deep sup
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename, '.session.mat']));

    deepSuperficialfromRipple = classification_DeepSuperficial(session);
end
