% custom_for_kenji

force_rerun = false;

re_run_if_before = '17-Aug-2021 11:00:00';

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

% if isempty(gcp('nocreate'))
%     parpool(4)
% end
WaitMessage = parfor_wait(length(sessions),'Waitbar',true);
% loop through each session
for i = 1:length(sessions)
    basename = sessions{i};
    basepath = [data_path,basename];
    disp(basepath)
    
    if ~isempty(re_run_if_before) &&...
            exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file')
        pass = check_file_date(fullfile(basepath,[basename,'.ripples.events.mat']),...
            re_run_if_before);
        if pass == 0
            WaitMessage.Send;
            continue
        end
    else
        pass=0;
    end
    
    % check for needed files
    if exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file') &&...
            ~force_rerun &&...
            ~pass ||...
            ~exist(basepath,'dir')
        WaitMessage.Send;
    else
        run_all(basepath,basename,force_rerun,pass)
        WaitMessage.Send;
    end
end
WaitMessage.Destroy;

custom_for_kenji_SWRunitMetrics()

function run_all(basepath,basename,force_rerun,pass)
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
if ~exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file') ||...
        force_rerun ||...
        pass
    addpath(genpath('D:\github\buzcode')) % cell-explorer and buzcode don't like eachother
    run_ripple_pipe_kenji(basepath,basename,spikes)
end

% make cell metrics
if ~exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') ||...
        force_rerun ||...
        pass
    rmpath(genpath('D:\github\buzcode')) % cell-explorer and buzcode don't like eachother
    
    cell_metrics = ProcessCellMetrics('basepath',basepath,...
        'showGUI',false,...
        'spikes',spikes,...
        'getWaveformsFromDat',false,...
        'manualAdjustMonoSyn',false,...
        'session',session);
    
    addpath(genpath('D:\github\buzcode'))
    
    close all
else
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
end

% assign region from kenji metadata
if contains(basepath,'Kenji') &&...
        ~isfield(cell_metrics,'brainRegion') ||...
        sum(contains(cell_metrics.brainRegion,'Unknown')) == length(cell_metrics.brainRegion) ||...
        force_rerun || pass
    try
        cell_metrics.brainRegion = get_kenji_region('basepath',basepath,...
            'basename',basename);
    catch
        warning([basepath,' get_kenji_region unit count ISSUE'])
        cell_metrics.brainRegion = get_kenji_region('basepath',basepath,...
            'basename',basename,'check_cell_count',false);
    end
    save(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics')
end
end

function pass = check_file_date(file,date_)
% check if file was ran before or on date
% true if file creation date is less than date_

files = dir(file);
[Y1, M1, D1, H1, ~, ~] = datevec(files.datenum);
[Y2, M2, D2, H2, ~, ~] = datevec(datenum(date_));
if D1==D2
    pass = H1 < H2;
else
    pass = D1 < D2;
end
end
% function run_ripple_pipe(basepath,basename,spikes)
%
% % load electrode position info
% load('A:\Data\Kenji\ElePosition.mat')
% shank_region = ElePosition(contains(ElePosition(:,2),basename),6:end);
% for i = 1:length(shank_region);shank_region{i}=lower(shank_region{i});end
%
% load(fullfile(basepath,[basename,'.session.mat']))
%
% % check if emg exist
% if ~exist(fullfile(basepath,[basename, '.EMGFromLFP.LFP.mat']),'file')
%     bz_EMGFromLFP(basepath,'basename',basename,...
%         'samplingFrequency',10,'savemat',true,'noPrompts',true);
% end
% % load lfp
% lfp = bz_GetLFP('all','basepath',basepath,'basename',basename,'noPrompts',true);
%
% % if data has ca1 then we can look at ripples that co-occur with sharp waves
% if any(strcmp(shank_region,'ca1') | strcmp(shank_region,'ca1c'))
%
%     % find out what channel has ca1 and run that channel
%     idx = find(strcmp(shank_region,'ca1') | strcmp(shank_region,'ca1c'));
%
%     chan_to_check = [session.extracellular.electrodeGroups.channels{idx}];
%     cell_idx = ismember(spikes.shankID,idx);
%
%     parfor i = 1:length(chan_to_check)
%         ripples = detect_ripple(lfp.data(:,chan_to_check(i)),...
%             lfp.timestamps,...
%             restrict_spikes(spikes,cell_idx));
%         n_rips(i) = length(ripples.peaks);
%     end
%     [~,c_idx] = max(n_rips);
%     c = chan_to_check(c_idx);
%
%     ripples = detect_ripple(lfp.data(:,c),...
%                             lfp.timestamps,...
%                             restrict_spikes(spikes,cell_idx));
%
%     ripples.detectorinfo.detectionparms.ripple_channel = c;
%
% elseif any(strcmp(shank_region,'ca3') |...
%         strcmp(shank_region,'ca2') |...
%         strcmp(shank_region,'ca') |...
%         strcmp(shank_region,'dg') |...
%         strcmp(shank_region,'dgca3'))
%
%     % find out what channel has HPC and run that channel
%     idx = find(strcmp(shank_region,'ca3') |...
%         strcmp(shank_region,'ca2') |...
%         strcmp(shank_region,'ca') |...
%         strcmp(shank_region,'dg') |...
%         strcmp(shank_region,'dgca3'));
%
%     chan_to_check = [session.extracellular.electrodeGroups.channels{idx}];
%     cell_idx = ismember(spikes.shankID,idx);
%
%     parfor i = 1:length(chan_to_check)
%
%         ripples = detect_ripple(lfp.data(:,chan_to_check(i)),...
%             lfp.timestamps,...
%             restrict_spikes(spikes,cell_idx));
%
%         n_rips(i) = length(ripples.peaks);
%     end
%     [~,c_idx] = max(n_rips);
%     c = chan_to_check(c_idx);
%
%     ripples = detect_ripple(lfp.data(:,c),...
%                             lfp.timestamps,...
%                             restrict_spikes(spikes,cell_idx));
%
%     ripples.detectorinfo.detectionparms.ripple_channel = c;
% else
%     disp('no HPC channels... find a channel you want to use and manually run')
%     keyboard
% end
%
% save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
% end
%
% function ripples = detect_ripple(lfp,timestamps,spikes)
%
% ripples = bz_FindRipples(lfp,timestamps,...
%                         'durations',[50,300],...
%                         'thresholds',[1,3],...
%                         'passband',[100,250],...
%                         'minDuration',20);
%
% ripples = eventSpikingTreshold(ripples,...
%                             'spikes',spikes,...
%                             'spikingThreshold',0.5,...
%                             'figOpt',false);
% end
% function temp_spikes = restrict_spikes(spikes,cell_idx)
% temp_spikes = spikes;
% n_cells = temp_spikes.numcells;
%
% for f = fields(temp_spikes)'
%     f = f{:};
%     if [isnumeric(temp_spikes.(f)) ||...
%             iscell(temp_spikes.(f))] &&...
%             length(temp_spikes.(f)) == n_cells
%         temp_spikes.(f) = temp_spikes.(f)(cell_idx);
%     end
% end
% end