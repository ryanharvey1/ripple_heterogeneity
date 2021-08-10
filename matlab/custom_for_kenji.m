% custom_for_kenji

load('A:\Data\Kenji\ElePosition.mat')
shank_region = ElePosition(:,6:end);
for i = 1:size(shank_region,1)
    for j = 1:size(shank_region,2)  
        shank_region{i,j}=lower(shank_region{i,j});
    end
end

idx = any(strcmp(shank_region,'ca1') |...
    strcmp(shank_region,'ca1c') |...
    strcmp(shank_region,'ca') |...
    strcmp(shank_region,'ca3') |...
    strcmp(shank_region,'ca2') |...
    strcmp(shank_region,'dg') |...
    strcmp(shank_region,'dgca3'),2);
                
sessions = ElePosition(idx,2);
    
% sessions ={...
%     'ec013.895_902','ec013.906_918','ec013.921_927',...
%     'ec013.931_942','ec013.944_958','ec013.961_974',...
%     'ec013.976_985','ec016.659_674','ec016.682_688',...
%     'ec016.694_711','ec016.715_735','ec016.740_764',...
%     'ec016.769_789','ec016.791_810','ec016.813_831',...
%     'ec016.835_850','ec016.853_867','ec016.871_889',...
%     'ec016.893_911','ec016.914_932','ec016.934_946',...
%     'ec016.950_965','ec016.969_986','ec016.1002_1023',...
%     'ec016.1025_1048','2006-4-18','2006-4-10',...
%     '2006-6-12','2006-6-13'...
%     };

data_path = 'A:\Data\Kenji\';

% loop through each session
for i = 1:length(sessions)
    basename = sessions{i};
    basepath = [data_path,basename];
    disp(basepath)
    % check for needed files
    if exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.SWRunitMetrics.mat']),'file')
        continue
    end
    
    % check and make session.mat
    if ~exist([basename '.session.mat'],'file')
        session = sessionTemplate(basepath,'showGUI',false);
        session.epochs = get_kenji_epochs('basepath',basepath,'basename',basename);
        save(fullfile(basepath,[basename '.session.mat']),'session')
    else
        load(fullfile(basepath,[basename '.session.mat']))
    end
    
    % load and process spikes
    spikes = loadSpikes(...
        'basepath',basepath,...
        'basename',basename,...
        'saveMat',true,...
        'format','Klustakwik',...
        'getWaveformsFromDat',false,...
        'getWaveformsFromSource',true,...
        'forceReload',false);
    
    % detect ripples
    if ~exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file')
        addpath(genpath('D:\github\buzcode')) % cell-explorer and buzcode don't like eachother
        run_ripple_pipe(basepath,basename,spikes)
    end
    
    % make cell metrics
    if ~exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file')
        rmpath(genpath('D:\github\buzcode')) % cell-explorer and buzcode don't like eachother
        
        cell_metrics = ProcessCellMetrics('basepath',basepath,...
            'showGUI',false,...
            'spikes',spikes,...
            'getWaveformsFromDat',false);
        
        addpath(genpath('D:\github\buzcode'))
        
        close all
    else
        load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    end
    
    % assign region from kenji metadata
    if contains(basepath,'Kenji') &&...
            ~isfield(cell_metrics,'brainRegion') ||...
            sum(contains(cell_metrics.brainRegion,'Unknown')) == length(cell_metrics.brainRegion)
        
        cell_metrics.brainRegion = get_kenji_region('basepath',basepath,...
            'basename',basename);
        save(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics')
    end
    
    % parse the epochs in order to make consistent pre/task/post structure
    % for swr code
    load(fullfile(basepath,[basename,'.ripples.events.mat']))
    
    % keep this at top of path
    addpath('D:\github\ripple_heterogeneity\matlab')
    parse_pre_task_post(session,basepath,basename,ripples,spikes)
end

function run_ripple_pipe(basepath,basename,spikes)

% regions = get_kenji_region('basepath',basepath,...
%     'basename',basename);

% load electrode position info
load('A:\Data\Kenji\ElePosition.mat')
shank_region = ElePosition(contains(ElePosition(:,2),basename),6:end);
for i = 1:length(shank_region);shank_region{i}=lower(shank_region{i});end

% if data has ca1 then we can look at ripples that co-occur with sharp waves
if any(strcmp(shank_region,'ca1') |...
        strcmp(shank_region,'ca') |...
        strcmp(shank_region,'ca1c'))
    disp('finding best ripple and sharp wave channels')
    if ~exist(fullfile(basepath,[basename '.swrCh.mat']),'file')
        swrCh = bz_swrChannels('basepath',basepath,'basename',basename,'Manual',true);
        save(fullfile(basepath,[basename '.swrCh.mat']),'swrCh')
    else
        load(fullfile(basepath,[basename '.swrCh.mat']))
    end
    
    if ~exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
        disp('detecting ripples')
        ripples = bz_DetectSWR([swrCh.ripple swrCh.sharpwave],'saveMat',true,'basename',basename);
        save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
    else
        load(fullfile(basepath,[basename '.ripples.events.mat']))
    end
elseif any(strcmp(shank_region,'ca3') |...
        strcmp(shank_region,'ca2') |...
        strcmp(shank_region,'dg') |...
        strcmp(shank_region,'dgca3'))

    % find out what channel has ca3 and run that channel
    idx = find(strcmp(shank_region,'ca3') |...
        strcmp(shank_region,'ca2') |...
        strcmp(shank_region,'dg') |...
        strcmp(shank_region,'dgca3'));
    
    load(fullfile(basepath,[basename '.sessionInfo.mat']))
    
    % test which channel within the spike group has the most detected ripples
    lfp = bz_GetLFP('all','basepath',basepath,'basename',basename,'noPrompts',true);
    i = 1;
    for c=sessionInfo.spikeGroups.groups{idx(1)}
        [ripples] = bz_FindRipples(lfp.data(:,c+1),lfp.timestamps);
        ripples = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.5);
        n_rips(i) = length(ripples.peaks);
        i=i+1;
    end
    [~,c_idx] = max(n_rips);
    c = sessionInfo.spikeGroups.groups{idx(1)}(c_idx);
    [ripples] = bz_FindRipples(lfp.data(:,c+1),lfp.timestamps);
        
    save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
else
%     lfp = bz_GetLFP('all','basepath',basepath,'basename',basename);
%     c = bz_GetBestRippleChan(lfp);
    disp('not sure how this session got through')
    c = input('pick channel number for ripples (0 indexing): ');
    [ripples] = bz_FindRipples(lfp.data(:,c+1),lfp.timestamps);
    save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end

load(fullfile(basepath,[basename,'.ripples.events.mat']))

disp('refining ripples by mua')
ripples = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.5);
save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end

