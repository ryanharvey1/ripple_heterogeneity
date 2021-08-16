% custom_for_kenji

force_rerun = false;

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
            ~force_rerun
    else
        run_all(basepath,basename,force_rerun)
    end
end

custom_for_kenji_SWRunitMetrics()

function run_all(basepath,basename,force_rerun)    
    % check and make session.mat
    if ~exist([basename '.session.mat'],'file')
        session = sessionTemplate(basepath,'showGUI',false);
        session.epochs = get_kenji_epochs('basepath',basepath,...
            'basename',basename);
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
    if ~exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file') || force_rerun
        addpath(genpath('D:\github\buzcode')) % cell-explorer and buzcode don't like eachother
        run_ripple_pipe(basepath,basename,spikes)
    end
    
    % make cell metrics
    if ~exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') || force_rerun
        rmpath(genpath('D:\github\buzcode')) % cell-explorer and buzcode don't like eachother
        
        cell_metrics = ProcessCellMetrics('basepath',basepath,...
            'showGUI',false,...
            'spikes',spikes,...
            'getWaveformsFromDat',false,...
            'manualAdjustMonoSyn',false,...
            'session',session,...
            'excludeMetrics',{'deepSuperficial'});
        
        addpath(genpath('D:\github\buzcode'))
        
        close all
    else
        load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    end
    
    % assign region from kenji metadata
    if contains(basepath,'Kenji') &&...
            ~isfield(cell_metrics,'brainRegion') ||...
            sum(contains(cell_metrics.brainRegion,'Unknown')) == length(cell_metrics.brainRegion) ||...
            force_rerun
        
        cell_metrics.brainRegion = get_kenji_region('basepath',basepath,...
                                                    'basename',basename);
        save(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics')
    end
    
    % parse the epochs in order to make consistent pre/task/post structure
    % for swr code
%     load(fullfile(basepath,[basename,'.ripples.events.mat']))
    
    % keep this at top of path
%     addpath('D:\github\ripple_heterogeneity\matlab')
%     parse_pre_task_post(session,basepath,basename,ripples,spikes)
end


% % keep this at top of path
% addpath('D:\github\ripple_heterogeneity\matlab')
% for i = 1:length(sessions)
%     basename = sessions{i};
%     basepath = [data_path,basename];
%     disp(basepath)
% 
%     if exist(fullfile(basepath,[basename,'.SWRunitMetrics.mat']),'file')
%         continue
%     end
%     
%     load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
% 
%     load(fullfile(basepath,[basename,'.ripples.events.mat']))
% 
%     load(fullfile(basepath,[basename '.session.mat']))
% 
%     parse_pre_task_post(session,basepath,basename,ripples,spikes)
%     
% end


function run_ripple_pipe(basepath,basename,spikes)

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
% load lfp
lfp = bz_GetLFP('all','basepath',basepath,'basename',basename,'noPrompts',true);

% if data has ca1 then we can look at ripples that co-occur with sharp waves
if any(strcmp(shank_region,'ca1') | strcmp(shank_region,'ca1c'))
    
    % find out what channel has ca1 and run that channel
    idx = find(strcmp(shank_region,'ca1') | strcmp(shank_region,'ca1c'));

    chan_to_check = [session.extracellular.electrodeGroups.channels{idx}];
    cell_idx = ismember(spikes.shankID,idx);

    parfor i = 1:length(chan_to_check)
        ripples = detect_ripple(lfp.data(:,chan_to_check(i)),...
            lfp.timestamps,...
            restrict_spikes(spikes,cell_idx));
        n_rips(i) = length(ripples.peaks);
    end
    [~,c_idx] = max(n_rips);
    c = chan_to_check(c_idx);
    
    ripples = detect_ripple(lfp.data(:,c),...
                            lfp.timestamps,...
                            restrict_spikes(spikes,cell_idx));
                        
    ripples.detectorinfo.detectionparms.ripple_channel = c;
    
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
        
    chan_to_check = [session.extracellular.electrodeGroups.channels{idx}];
    cell_idx = ismember(spikes.shankID,idx);

    parfor i = 1:length(chan_to_check)
        
        ripples = detect_ripple(lfp.data(:,chan_to_check(i)),...
            lfp.timestamps,...
            restrict_spikes(spikes,cell_idx));

        n_rips(i) = length(ripples.peaks);
    end
    [~,c_idx] = max(n_rips);
    c = chan_to_check(c_idx);
    
    ripples = detect_ripple(lfp.data(:,c),...
                            lfp.timestamps,...
                            restrict_spikes(spikes,cell_idx));
                        
    ripples.detectorinfo.detectionparms.ripple_channel = c;
else
    disp('no HPC channels... find a channel you want to use and manually run')
    keyboard
end

save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end

function ripples = detect_ripple(lfp,timestamps,spikes)

ripples = bz_FindRipples(lfp,timestamps,...
                        'durations',[50,300],...
                        'thresholds',[1,3],...
                        'passband',[100,250],...
                        'minDuration',20);
    
ripples = eventSpikingTreshold(ripples,...
                            'spikes',spikes,...
                            'spikingThreshold',0.5,...
                            'figOpt',false);
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