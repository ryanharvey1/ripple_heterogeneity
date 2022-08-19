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
% basepaths = {
%     'Z:\Data\FujisawaS\EE\EE0711',...
%     'Z:\Data\FujisawaS\GG\GG0401',...
%     'Z:\Data\FujisawaS\GG\GG0406'}
for i = 1:length(basepaths)
    basepath = basepaths{i};

    spikes = loadSpikes('basepath',basepath,...
        'getWaveformsFromSource',true,...
        'getWaveformsFromDat',true,...
        'format','klustakwik','forceReload',true);
end

for i = 1:length(basepaths)
    
end

%% get cell metrics
for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,...
        'manualAdjustMonoSyn',false,'showGUI',false,...
        'excludeManipulationIntervals',false,'forceReload',true);
end


%% detect ripples !!! Manual process !!!

basepaths = {'Z:\Data\FujisawaS\FF\FF1114',...
    'Z:\Data\FujisawaS\FF\FF1116',...
    'Z:\Data\FujisawaS\FF\FF1119',...
    'Z:\Data\FujisawaS\EE\EE0710',...
    'Z:\Data\FujisawaS\EE\EE0711',...
    'Z:\Data\FujisawaS\GG\GG0401',...
    'Z:\Data\FujisawaS\GG\GG0406'}
for i = 1:length(basepaths)
   cd(basepaths{i}) 
   NeuroScope2()
end
pause;
basepath = pwd
basename = basenameFromBasepath(basepath);

load(fullfile(basepath,[basename,'.session.mat']))

ripples = DetectSWR([session.channelTags.Ripple.channels,...
    session.channelTags.SharpWave.channels,...
    session.brainRegions.PFC.channels(1)],...
    'basepath',basepath,...
    'saveMat',true,'thresSDswD', [0.1, .5],'thresSDrip', [0.5, 1.5],...
    'forceDetect',true,'check',false);

ripples = FindRipples('basepath',basepath,...
    'channel',session.channelTags.Ripple.channels,'noise',16,'duration',[50,200],'saveMat',true);


% refine ripple detection using spiking level
spikes = importSpikes('basepath',basepath,'cellType', "Pyramidal Cell", 'brainRegion', "CA1");
ripplesTemp = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.4);
ripples = ripplesTemp;

save(fullfile(basepath,[basename,'.ripples.events.mat']),'ripples')

%% deep sup

basepaths = {'Z:\Data\FujisawaS\FF\FF1114',...
    'Z:\Data\FujisawaS\FF\FF1116',...
    'Z:\Data\FujisawaS\FF\FF1119',...
    'Z:\Data\FujisawaS\EE\EE0710',...
    'Z:\Data\FujisawaS\EE\EE0711',...
    'Z:\Data\FujisawaS\GG\GG0401',...
    'Z:\Data\FujisawaS\GG\GG0406'}
for i = 1:length(basepaths)
   cd(basepaths{i}) 
   gui_session()
end

for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename, '.session.mat']));

    deepSuperficialfromRipple = classification_DeepSuperficial(session);
end

%% update epochs
basepaths = {'Z:\Data\FujisawaS\FF\FF1114',...
    'Z:\Data\FujisawaS\FF\FF1116',...
    'Z:\Data\FujisawaS\FF\FF1119',...
    'Z:\Data\FujisawaS\EE\EE0710',...
    'Z:\Data\FujisawaS\EE\EE0711',...
    'Z:\Data\FujisawaS\GG\GG0401',...
    'Z:\Data\FujisawaS\GG\GG0406'}
for i = 1:length(basepaths)
   cd(basepaths{i}) 
   gui_session()
end

%% add ripple metrics to cell_metrics

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

    cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,...
        'manualAdjustMonoSyn',false,'showGUI',false,'getWaveformsFromDat',false);
end

%%
% pull out basepaths and basenames
basename = []
for i = 1:length(basepaths)
    basename{i} = basenameFromBasepath(basepaths{i});
end
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basename);
cell_metrics = CellExplorer('metrics',cell_metrics);

%%

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    general_behavior_file('basepath',basepath)
end